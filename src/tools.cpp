// lq export / matrix -- outside the quantification pipeline
//
//   export   binary outputs -> the text files of LineageQuant 1.0.0, byte for byte
//              lq export S.lqres [-o DIR] [--name S] [--em-name EM_Final_Expression.tsv]
//                -> DIR/EM_Final_Expression.tsv, DIR/S.idxstats.tsv, DIR/S.efflen.tsv
//              lq export Global_Theta_Matrix.lqb [-o FILE.tsv]
//   matrix   genome-wide expression matrix (port of core/aggregator.py, Step 04 of 1.0.0; same output)
//              lq matrix -r RESULTS_DIR -o ORPHANS [--out LineageQuant_Final_Merged_TPM.tsv]
//                                  [--mode allo|auto] [--tpm-length full|eff]
//            reads <sample>.lqres (or the text files of 1.0.0: <sample>.idxstats.tsv + EM_Final_*.tsv
//            [+ <sample>.efflen.tsv]) found anywhere below RESULTS_DIR.
#include "binio.hpp"
#include "common.hpp"

#include <dirent.h>
#include <functional>
#include <map>
#include <set>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>

namespace lq {

namespace {

std::string base_name(const std::string& p) {
    size_t s = p.find_last_of('/');
    return s == std::string::npos ? p : p.substr(s + 1);
}
std::string dir_name(const std::string& p) {
    size_t s = p.find_last_of('/');
    if (s == std::string::npos) return ".";
    return s == 0 ? "/" : p.substr(0, s);
}
std::string strip_suffix(const std::string& s, const std::string& x) {
    return ends_with(s, x) ? s.substr(0, s.size() - x.size()) : s;
}

void export_usage() {
    std::fprintf(stderr,
                 "usage: lq export S.lqres [-o DIR] [--name S] [--em-name EM_Final_Expression.tsv]\n"
                 "       lq export Global_Theta_Matrix.lqb [-o FILE]\n"
                 "  writes the text files of LineageQuant 1.0.0 (same bytes as the text writers):\n"
                 "    .lqres -> DIR/<em-name>, DIR/S.idxstats.tsv, DIR/S.efflen.tsv  (DIR default: next to the file)\n"
                 "    .lqb   -> Global_Theta_Matrix.tsv (default: next to the file)\n");
}

}  // namespace

int export_main(int argc, char** argv) {
    std::string in, out, name, em_name = "EM_Final_Expression.tsv";
    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&]() -> std::string {
            if (i + 1 >= argc) { export_usage(); die("missing value for " + a); }
            return argv[++i];
        };
        if (a == "-o" || a == "--out") out = need();
        else if (a == "--name") name = need();
        else if (a == "--em-name") em_name = need();
        else if (a == "-h" || a == "--help") { export_usage(); return 0; }
        else if (!a.empty() && a[0] == '-') { export_usage(); die("unknown option " + a); }
        else if (in.empty()) in = a;
        else { export_usage(); die("one input file at a time"); }
    }
    if (in.empty()) { export_usage(); return 1; }
    if (has_magic(in, THETA_MAGIC)) {
        ThetaBin T;
        read_theta_bin(in, T);
        if (out.empty()) out = dir_name(in) + "/Global_Theta_Matrix.tsv";
        write_text_file(out, [&](FILE* f) {
            theta_tsv_header(f);
            for (auto& r : T.rows) theta_tsv_row(f, T, r);
        });
        logm("[SUCCESS] " + out + " (" + commas((long long)T.rows.size()) + " copies)");
        return 0;
    }
    if (!has_magic(in, RES_MAGIC)) die(in + ": neither a .lqres result nor a .lqb Theta file");
    ResultBin R;
    read_result_bin(in, R);
    if (name.empty()) name = strip_suffix(base_name(in), ".lqres");
    if (out.empty()) out = dir_name(in);
    mkdirs(out);
    const std::string em = out + "/" + em_name, idx = out + "/" + name + ".idxstats.tsv", eff = out + "/" + name + ".efflen.tsv";
    write_text_file(em, [&](FILE* f) { write_em_tsv(f, R); });
    write_text_file(idx, [&](FILE* f) { write_idx_tsv(f, R); });
    write_text_file(eff, [&](FILE* f) { write_eff_tsv(f, R); });
    logm("[SUCCESS] " + em + ", " + idx + ", " + eff + " (" + commas((long long)R.T()) + " targets, " +
         commas((long long)R.refs.size()) + " references)");
    return 0;
}

// ===================================================================================== matrix
namespace {

void matrix_usage() {
    std::fprintf(stderr,
                 "usage: lq matrix -r RESULTS_DIR -o ORPHANS [--out LineageQuant_Final_Merged_TPM.tsv]\n"
                 "                           [--mode allo|auto|idxstats] [--tpm-length full|eff]\n"
                 "  genome-wide matrix: per sample <sample>_RawReads (EM-assigned fragments) and <sample>_TPM;\n"
                 "  same file as core/aggregator.py of LineageQuant 1.0.0. Not part of `lq run`.\n"
                 "  --tpm-length full (default): reads per kb of the longest transcript; eff: reads per effective\n"
                 "  length used by the EM (genes without one fall back to full length).\n"
                 "  --mode idxstats (ablation, as 1.0.0): every gene from the per-reference counts, i.e. every reported\n"
                 "  alignment counts in full (featureCounts -M -O); none of LineageQuant's algorithms.\n"
                 "  --out is relative to RESULTS_DIR unless absolute.\n");
}

// Python str.strip() for str.isspace() ASCII characters
std::string_view pstrip(std::string_view s) {
    auto ws = [](unsigned char c) { return c == ' ' || (c >= 9 && c <= 13) || (c >= 0x1c && c <= 0x1f); };
    size_t a = 0, b = s.size();
    while (a < b && ws((unsigned char)s[a])) ++a;
    while (b > a && ws((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}

// glob.glob(os.path.join(root, "**", pattern), recursive=True) order: files of a directory in readdir order,
// then each subdirectory (pre-order, readdir order); names starting with "." are skipped; symlinks followed
void walk(const std::string& dir, int depth, std::set<std::pair<dev_t, ino_t>>& seen,
          const std::function<void(const std::string&, const std::string&, int)>& visit) {
    struct stat st;
    if (stat(dir.c_str(), &st) != 0 || !S_ISDIR(st.st_mode)) return;
    if (!seen.insert({st.st_dev, st.st_ino}).second) return;
    DIR* d = opendir(dir.c_str());
    if (!d) return;
    std::vector<std::pair<std::string, bool>> ents;
    while (dirent* e = readdir(d)) {
        std::string n = e->d_name;
        if (n.empty() || n[0] == '.') continue;
        std::string p = dir + "/" + n;
        struct stat s2;
        bool isdir = stat(p.c_str(), &s2) == 0 && S_ISDIR(s2.st_mode);
        ents.emplace_back(n, isdir);
    }
    closedir(d);
    for (auto& e : ents) if (!e.second) visit(dir, e.first, depth);
    for (auto& e : ents) if (e.second) walk(dir + "/" + e.first, depth + 1, seen, visit);
}

struct SampleSrc {
    std::string path;     // .lqres or .idxstats.tsv
    bool bin = false;
    int depth = 0;
};

// pandas read_csv(sep="\t") of an EM_Final_*.tsv: the three columns the aggregator uses
void read_em_tsv(const std::string& path, std::vector<std::string>& allele, std::vector<std::string>& cluster,
                 std::vector<double>& reads) {
    LineReader rd(path);
    std::string_view l;
    if (!rd.next(l)) die(path + ": empty EM file");
    auto h = split_tab(l);
    int ia = -1, ic = -1, ir = -1;
    for (size_t i = 0; i < h.size(); ++i) {
        if (h[i] == "AlleleID") ia = (int)i;
        else if (h[i] == "ClusterID") ic = (int)i;
        else if (h[i] == "EM_Allocated_Reads") ir = (int)i;
    }
    if (ia < 0 || ic < 0 || ir < 0) die(path + ": needs columns ClusterID, AlleleID, EM_Allocated_Reads");
    while (rd.next(l)) {
        if (l.empty()) continue;   // pandas skips blank lines
        auto p = split_tab(l);
        auto get = [&](int i) { return i < (int)p.size() ? std::string(p[i]) : std::string(); };
        allele.push_back(get(ia));
        cluster.push_back(get(ic));
        double v = 0.0;
        if (!py_float(get(ir), v)) v = std::nan("");
        reads.push_back(v);
    }
}

}  // namespace

int matrix_main(int argc, char** argv) {
    std::string results, orphan, out = "LineageQuant_Final_Merged_TPM.tsv", mode = "allo", tpm_length = "full";
    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&]() -> std::string {
            if (i + 1 >= argc) { matrix_usage(); die("missing value for " + a); }
            return argv[++i];
        };
        if (a == "-r" || a == "--results_dir" || a == "--results-dir") results = need();
        else if (a == "-o" || a == "--orphan") orphan = need();
        else if (a == "--out") out = need();
        else if (a == "--mode") mode = need();
        else if (a == "--tpm-length" || a == "--tpm_length") tpm_length = need();
        else if (a == "-h" || a == "--help") { matrix_usage(); return 0; }
        else { matrix_usage(); die("unknown argument " + a); }
    }
    if (results.empty() || orphan.empty()) { matrix_usage(); die("-r and -o are required"); }
    if (mode != "allo" && mode != "auto" && mode != "idxstats") die("--mode must be allo, auto or idxstats");
    const bool use_em = mode != "idxstats";
    if (tpm_length != "full" && tpm_length != "eff") die("--tpm-length must be full or eff");
    const bool by_eff = tpm_length == "eff" && mode != "idxstats";
    const double T0 = wall();
    while (results.size() > 1 && results.back() == '/') results.pop_back();

    // ---- orphan genes
    std::unordered_set<std::string> orphan_set;
    {
        LineReader rd(orphan);
        std::string_view l;
        while (rd.next(l)) {
            std::string_view s = pstrip(l);
            if (s.empty() || (!l.empty() && l[0] == '#')) continue;
            auto parts = split_tab(s);
            if (parts.size() >= 2) {
                std::string_view g = parts[1];
                size_t b = 0;
                for (;;) {
                    size_t e = g.find(',', b);
                    std::string_view x = g.substr(b, e == std::string_view::npos ? std::string_view::npos : e - b);
                    orphan_set.insert(std::string(rsplit_dot(x)));
                    if (e == std::string_view::npos) break;
                    b = e + 1;
                }
            } else orphan_set.insert(std::string(rsplit_dot(parts[0])));
        }
    }
    logm("[INFO] orphan genes: " + commas((long long)orphan_set.size()) + " (" + orphan + ")");

    // ---- samples: <sample>.lqres or <sample>.idxstats.tsv, the shallowest one per sample name
    std::map<std::string, SampleSrc> src;
    std::set<std::pair<dev_t, ino_t>> seen;
    walk(results, 0, seen, [&](const std::string& dir, const std::string& n, int depth) {
        bool bin = ends_with(n, ".lqres"), txt = ends_with(n, ".idxstats.tsv");
        if (!bin && !txt) return;
        std::string s = bin ? strip_suffix(n, ".lqres") : strip_suffix(n, ".idxstats.tsv");
        auto it = src.find(s);
        if (it == src.end() || depth < it->second.depth || (depth == it->second.depth && bin && !it->second.bin))
            src[s] = SampleSrc{dir + "/" + n, bin, depth};
    });
    if (src.empty()) {
        logm("[ERROR] no <sample>.lqres (lq quant) or <sample>.idxstats.tsv found below " + results);
        return 1;
    }
    std::vector<std::string> samples;
    for (auto& kv : src) samples.push_back(kv.first);   // std::map: sorted like Python sorted()
    logm("[INFO] " + std::to_string(samples.size()) + " samples. Mode=" + mode + " (EM=yes). Aggregating...");

    // ---- per gene, per sample (insertion order = the aggregator's dict order: it decides the TPM sums)
    std::unordered_map<std::string, uint32_t> gix;
    std::vector<std::string> genes;
    std::vector<std::vector<std::pair<uint32_t, double>>> counts;   // gene -> (sample, reads), last write wins
    std::vector<std::vector<std::pair<uint32_t, double>>> rate;     // --tpm-length eff
    std::vector<long long> glen;                                     // -1 = none
    std::vector<int> gsrc;                                           // -1 none, 0 idxstats, 1 em, 2 em_orphan
    auto gene_id = [&](const std::string& g) -> uint32_t {
        auto it = gix.find(g);
        if (it != gix.end()) return it->second;
        uint32_t id = (uint32_t)genes.size();
        gix.emplace(g, id);
        genes.push_back(g);
        counts.emplace_back();
        rate.emplace_back();
        glen.push_back(-1);
        gsrc.push_back(-1);
        return id;
    };
    auto set_count = [&](uint32_t g, uint32_t s, double v) {
        for (auto& x : counts[g]) if (x.first == s) { x.second = v; return; }
        counts[g].emplace_back(s, v);
    };
    // gene ids seen in the length table (gene_lengths), without creating rows
    std::unordered_map<std::string, long long> gene_lengths;
    size_t n_missing_em = 0, n_no_eff = 0;
    for (uint32_t si = 0; si < samples.size(); ++si) {
        const std::string& sample = samples[si];
        const SampleSrc& S = src[sample];
        ResultBin R;
        std::vector<std::string> allele, cluster;
        std::vector<double> reads;
        bool has_em = false;
        std::unordered_map<std::string, double> eff;
        // idxstats
        std::vector<std::string> ord;                       // idx_counts insertion order
        std::unordered_map<std::string, double> idx_counts;
        auto add_idx = [&](const std::string& ref, long long len, unsigned long long mapped) {
            std::string g(rsplit_dot(ref));
            auto it = gene_lengths.find(g);
            if (it == gene_lengths.end() || len > it->second) gene_lengths[g] = len;
            auto jt = idx_counts.find(g);
            if (jt == idx_counts.end()) { ord.push_back(g); idx_counts.emplace(g, 0.0 + (double)mapped / 2.0); }
            else jt->second = jt->second + (double)mapped / 2.0;
        };
        if (!use_em) {   // ablation: idxstats only
            if (S.bin) {
                read_result_bin(S.path, R);
                for (size_t i = 0; i < R.refs.size(); ++i)
                    if (R.refs[i] != "*") add_idx(R.refs[i], R.rlen[i], R.mapped[i]);
            } else {
                LineReader rd(S.path);
                std::string_view l;
                while (rd.next(l)) {
                    auto p = split_tab(l);
                    if (p.size() >= 4 && p[0] != "*") add_idx(std::string(p[0]), parse_ll(p[1]), (unsigned long long)parse_ll(p[2]));
                }
            }
            for (auto& g : ord) {
                uint32_t id = gene_id(g);
                set_count(id, si, idx_counts[g]);
                if (gsrc[id] < 0) gsrc[id] = 0;
            }
            continue;
        }
        if (S.bin) {
            read_result_bin(S.path, R);
            for (size_t i = 0; i < R.refs.size(); ++i)
                if (R.refs[i] != "*") add_idx(R.refs[i], R.rlen[i], R.mapped[i]);
            has_em = true;
            allele = R.allele;
            cluster = R.cluster;
            reads.resize(R.T());
            for (size_t i = 0; i < R.T(); ++i) reads[i] = py_round(R.counts[i], 2);      // what the text holds
            if (by_eff)
                for (size_t i = 0; i < R.T(); ++i) eff[R.allele[i]] = py_round(R.eff[i], 4);
        } else {
            LineReader rd(S.path);
            std::string_view l;
            while (rd.next(l)) {
                auto p = split_tab(l);
                if (p.size() >= 4 && p[0] != "*") add_idx(std::string(p[0]), parse_ll(p[1]), (unsigned long long)parse_ll(p[2]));
            }
            const std::string dir = dir_name(S.path);
            for (const char* nm : {"EM_Final_Hybrid.tsv", "EM_Final_Expression_v2.tsv", "EM_Final_Expression.tsv"}) {
                std::string p = dir + "/" + nm;
                if (file_exists(p)) { read_em_tsv(p, allele, cluster, reads); has_em = true; break; }
            }
            if (has_em && by_eff) {
                std::string p = dir + "/" + sample + ".efflen.tsv";
                if (file_exists(p)) {
                    LineReader er(p);
                    bool first = true;
                    while (er.next(l)) {
                        if (first) { first = false; continue; }
                        auto x = split_tab(l);
                        double v;
                        if (x.size() >= 3) { if (!py_float(x[2], v)) die(p + ": not a number: " + std::string(x[2])); eff[std::string(x[0])] = v; }
                    }
                }
            }
        }
        if (!has_em) {
            ++n_missing_em;
            logm("[WARN] " + sample + ": no EM output found; multi-copy genes will be missing for this sample (did quant finish?)");
            for (auto& g : ord) {
                if (!orphan_set.count(g)) continue;
                uint32_t id = gene_id(g);
                set_count(id, si, idx_counts[g]);
                if (gsrc[id] < 0) gsrc[id] = 0;
            }
            continue;
        }
        // groupby(Clean_Gene)["EM_Allocated_Reads"].sum(): Kahan per group, groups in sorted order
        std::map<std::string, std::pair<double, double>> em_sum;   // gene -> (sum, compensation)
        std::unordered_set<std::string> orph_em;
        std::vector<std::string> clean(allele.size());
        for (size_t i = 0; i < allele.size(); ++i) {
            clean[i] = std::string(rsplit_dot(allele[i]));
            auto& sc = em_sum[clean[i]];
            const double val = reads[i];
            if (std::isnan(val)) continue;
            double y = val - sc.second;
            double t = sc.first + y;
            sc.second = (t - sc.first) - y;
            if (sc.second != sc.second) sc.second = 0;
            sc.first = t;
        }
        for (size_t i = 0; i < allele.size(); ++i) if (cluster[i] == "ORPHAN") orph_em.insert(clean[i]);
        for (auto& kv : em_sum) {
            uint32_t id = gene_id(kv.first);
            set_count(id, si, kv.second.first);
            gsrc[id] = orph_em.count(kv.first) ? 2 : 1;
        }
        if (by_eff) {
            if (eff.empty()) {
                ++n_no_eff;
                logm("[WARN] " + sample + ": no " + sample + ".efflen.tsv -> TPM by full length for this sample");
            } else {
                for (size_t i = 0; i < allele.size(); ++i) {
                    auto it = eff.find(allele[i]);
                    if (it == eff.end() || !(it->second > 0)) continue;   // `if e and e > 0`
                    uint32_t id = gix[clean[i]];
                    auto& rv = rate[id];
                    bool done = false;
                    for (auto& x : rv) if (x.first == si) { x.second = x.second + reads[i] / it->second; done = true; break; }
                    if (!done) rv.emplace_back(si, 0.0 + reads[i] / it->second);
                }
            }
        }
        for (auto& g : ord) {
            if (em_sum.count(g)) continue;
            uint32_t id = gene_id(g);
            set_count(id, si, idx_counts[g]);
            if (gsrc[id] < 0) gsrc[id] = 0;
        }
    }
    if (use_em && n_missing_em == samples.size()) {
        logm("[ERROR] no sample has EM output (lq quant).");
        return 1;
    }
    (void)n_no_eff;

    // ---- RawReads / RPK / TPM
    logm(std::string("[INFO] Computing TPM normalization") + (by_eff ? " (effective lengths)" : "") + "...");
    const size_t G = genes.size(), NS = samples.size();
    std::vector<double> raw(G * NS, 0.0), rpk(NS * G, 0.0);
    std::vector<long long> len(G);
    for (size_t g = 0; g < G; ++g) {
        auto it = gene_lengths.find(genes[g]);
        len[g] = it == gene_lengths.end() ? 1000 : it->second;
        const double kb = (double)len[g] / 1000.0;
        for (size_t s = 0; s < NS; ++s) {
            double r = 0.0;
            for (auto& x : counts[g]) if (x.first == s) r = x.second;
            raw[g * NS + s] = r;
            const std::pair<uint32_t, double>* rt = nullptr;
            if (by_eff) for (auto& x : rate[g]) if (x.first == s) rt = &x;
            rpk[s * G + g] = rt ? rt->second * 1000.0 : (kb > 0 ? r / kb : 0.0);
        }
    }
    std::vector<double> tpm(NS * G);
    for (size_t s = 0; s < NS; ++s) {
        const double* col = rpk.data() + s * G;
        std::vector<double> nz(col, col + G);           // Series.sum(skipna=True): NaN -> 0
        for (auto& x : nz) if (std::isnan(x)) x = 0.0;
        const double total = np_add_reduce(nz.data(), G);
        const double scaling = total > 0 ? total / 1000000.0 : 1.0;
        for (size_t g = 0; g < G; ++g) tpm[s * G + g] = std::nearbyint((col[g] / scaling) * 100.0) / 100.0;
    }
    auto type_of = [&](size_t g) -> const char* {
        return gsrc[g] == 2 ? "Orphan(EM)" : gsrc[g] == 1 ? "Polyploid(EM)"
               : (use_em || orphan_set.count(genes[g])) ? "Orphan(BAM_Base)" : "Polyploid(idxstats)";
    };
    std::vector<uint32_t> order(G);
    for (size_t g = 0; g < G; ++g) order[g] = (uint32_t)g;
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
        int c = std::strcmp(type_of(a), type_of(b));
        if (c != 0) return c > 0;                 // Type descending
        return genes[a] < genes[b];               // GeneID ascending
    });
    // ---- write (pandas to_csv: sep="\t", index=False, floats as repr, csv minimal quoting)
    auto field = [](const std::string& s) {
        if (s.find_first_of("\t\"\n\r") == std::string::npos) return s;
        std::string o = "\"";
        for (char c : s) { if (c == '"') o += '"'; o += c; }
        return o + "\"";
    };
    auto num = [](double x) { return std::isnan(x) ? std::string() : py_repr(x); };   // na_rep=""
    const std::string out_path = (!out.empty() && out[0] == '/') ? out : results + "/" + out;
    size_t n_poly = 0;
    write_text_file(out_path, [&](FILE* f) {
        std::string h = "GeneID\tType\tLength_bp";
        for (auto& s : samples) h += "\t" + field(s + "_RawReads");
        for (auto& s : samples) h += "\t" + field(s + "_TPM");
        fputs(h.c_str(), f);
        fputc('\n', f);
        std::string line;
        for (uint32_t g : order) {
            line = field(genes[g]);
            line += '\t';
            line += type_of(g);
            line += '\t';
            line += std::to_string(len[g]);
            for (size_t s = 0; s < NS; ++s) { line += '\t'; line += num(raw[g * NS + s]); }
            for (size_t s = 0; s < NS; ++s) { line += '\t'; line += num(tpm[s * G + g]); }
            line += '\n';
            fwrite(line.data(), 1, line.size(), f);
            if (gsrc[g] == 1) ++n_poly;
        }
    });
    logm("[INFO] " + commas((long long)G) + " genes (" + commas((long long)n_poly) + " allele copies / " +
         commas((long long)(G - n_poly)) + " single-copy)");
    logm("[SUCCESS] TPM matrix: " + out_path + sfmt(" (%.1fs)", wall() - T0));
    return 0;
}

}  // namespace lq
