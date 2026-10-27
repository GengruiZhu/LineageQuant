// LineageQuant 1.1.0 -- reference preparation (C++ ports of the 1.0.0 Python modules; same outputs)
//
//   snpdb-msa            lineage_snp/build_snp_msa.py + export_snp_db.py --mode msa      (MAFFT)
//   snpdb-vcf            lineage_snp/build_snp_vcf.py + export_snp_db.py --mode population
//   build-markers        lineage_kmer/build_markers.py + filter_low_complexity.py + txt2fasta.py   (KMC)
//   extract-transcripts  utils/genome_pair.py                                             (gffread)
//   normalize-clusters   cluster/normalize_clusters.py
//
// The SNP databases are written directly in the text format lqcore reads (the pickle step of 1.0.0 is gone).
#include "prep.hpp"
#include "proc.hpp"

#include <functional>

#include <atomic>
#include <ftw.h>
#include <map>
#include <mutex>
#include <regex>
#include <set>
#include <thread>
#include <unordered_map>
#include <unordered_set>

namespace lq {

// ===================================================================================== small Python helpers
namespace {
inline bool py_space(unsigned char c) { return c == ' ' || (c >= 9 && c <= 13) || (c >= 0x1c && c <= 0x1f); }
}
std::string_view py_strip(std::string_view s) {
    size_t a = 0, b = s.size();
    while (a < b && py_space((unsigned char)s[a])) ++a;
    while (b > a && py_space((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}
std::string_view py_rstrip(std::string_view s) {
    size_t b = s.size();
    while (b > 0 && py_space((unsigned char)s[b - 1])) --b;
    return s.substr(0, b);
}
// str.split(sep) (sep non-empty)
std::vector<std::string_view> py_split(std::string_view s, char sep) {
    std::vector<std::string_view> out;
    size_t b = 0;
    for (;;) {
        size_t e = s.find(sep, b);
        if (e == std::string_view::npos) { out.emplace_back(s.substr(b)); return out; }
        out.emplace_back(s.substr(b, e - b));
        b = e + 1;
    }
}
std::vector<std::string_view> py_split_ws(std::string_view s) {
    std::vector<std::string_view> out;
    size_t i = 0, n = s.size();
    while (i < n) {
        while (i < n && py_space((unsigned char)s[i])) ++i;
        if (i >= n) break;
        size_t j = i;
        while (j < n && !py_space((unsigned char)s[j])) ++j;
        out.emplace_back(s.substr(i, j - i));
        i = j;
    }
    return out;
}
// s.split(".")[0]
std::string_view before_dot(std::string_view s) {
    size_t p = s.find('.');
    return p == std::string_view::npos ? s : s.substr(0, p);
}
// Python int(s): surrounding whitespace, optional sign, digits with single "_" between them
bool py_int(std::string_view s, long long& v) {
    s = py_strip(s);
    size_t i = 0;
    bool neg = false;
    if (i < s.size() && (s[i] == '+' || s[i] == '-')) neg = s[i++] == '-';
    if (i >= s.size()) return false;
    long long x = 0;
    bool prev_digit = false;
    for (; i < s.size(); ++i) {
        char c = s[i];
        if (c >= '0' && c <= '9') { x = x * 10 + (c - '0'); prev_digit = true; }
        else if (c == '_' && prev_digit && i + 1 < s.size() && isdigit((unsigned char)s[i + 1])) prev_digit = false;
        else return false;
    }
    v = neg ? -x : x;
    return true;
}
std::string lower_ascii(std::string_view s) {
    std::string o(s);
    for (auto& c : o) c = (char)tolower((unsigned char)c);
    return o;
}

// Bio.SeqIO.parse(path, "fasta"): id = first word of the title, sequence = lines right-stripped and joined,
// spaces and "\r" removed; text before the first ">" is ignored
void read_fasta_bio(const std::string& path, const std::function<void(std::string&&, std::string&&)>& rec) {
    LineReader rd(path);
    std::string_view l;
    bool in = false;
    std::string id, seq;
    auto flush = [&]() {
        std::string s;
        s.reserve(seq.size());
        for (char c : seq) if (c != ' ' && c != '\r') s += c;
        rec(std::move(id), std::move(s));
        id.clear();
        seq.clear();
    };
    while (rd.next(l)) {
        if (!l.empty() && l[0] == '>') {
            if (in) flush();
            in = true;
            auto w = py_split_ws(py_rstrip(l.substr(1)));
            id = w.empty() ? std::string() : std::string(w[0]);
            continue;
        }
        if (!in) continue;
        auto r = py_rstrip(l);
        seq.append(r.data(), r.size());
    }
    if (in) flush();
}

// ancestry file as build_snp_msa.py / utils.lineage read it: {gene.split(".")[0]: LINEAGE.upper()}
void read_ancestry_upper(const std::string& path, std::unordered_map<std::string, std::string>& anc, bool skip_hash) {
    LineReader rd(path);
    std::string_view l;
    while (rd.next(l)) {
        if (py_strip(l).empty() || l.rfind("GeneID", 0) == 0) continue;
        if (skip_hash && !py_strip(l).empty() && py_strip(l)[0] == '#') continue;
        auto p = py_split_ws(py_strip(l));
        if (p.size() >= 2) anc[std::string(before_dot(p[0]))] = upper(p[1]);
    }
}

static int rm_cb(const char* p, const struct stat*, int, struct FTW*) { return remove(p); }
void remove_dir_tree(const std::string& d) { nftw(d.c_str(), rm_cb, 64, FTW_DEPTH | FTW_PHYS); }

std::string make_temp_dir(const std::string& base, const std::string& prefix) {
    std::string root = base;
    if (root.empty()) { const char* t = getenv("TMPDIR"); root = t && *t ? t : "/tmp"; }
    mkdirs(root);
    std::string tmpl = root + "/" + prefix + "XXXXXX";
    std::vector<char> b(tmpl.begin(), tmpl.end());
    b.push_back('\0');
    if (!mkdtemp(b.data())) die("cannot create a temporary directory in " + root);
    return std::string(b.data());
}

// ===================================================================================== SNP DB writer (export_snp_db.write_lq)
namespace {
struct SnpDb {
    // insertion-ordered dict {(tid, pos): entry}; overwriting keeps the position, replaces the value
    std::vector<std::string> tids;
    std::unordered_map<std::string, uint32_t> tix;
    std::unordered_map<uint64_t, uint32_t> key;
    std::vector<uint32_t> e_tid, e_pos;
    std::vector<std::string> e_val;     // msa: expected base per lineage of the universe ('\0' = none); vcf: REF ALT AF...
    uint32_t tid_of(const std::string& t) {
        auto it = tix.find(t);
        if (it != tix.end()) return it->second;
        uint32_t i = (uint32_t)tids.size();
        tix.emplace(t, i);
        tids.push_back(t);
        return i;
    }
    void put(uint32_t tid, uint32_t pos, std::string&& v) {
        uint64_t k = ((uint64_t)tid << 32) | pos;
        auto it = key.find(k);
        if (it != key.end()) { e_val[it->second] = std::move(v); return; }
        key.emplace(k, (uint32_t)e_tid.size());
        e_tid.push_back(tid);
        e_pos.push_back(pos);
        e_val.push_back(std::move(v));
    }
    size_t size() const { return e_tid.size(); }
};
// grouped by transcript in the DB's order, positions ascending; row(i) prints the columns after tid and pos
void write_lq(const SnpDb& db, const std::string& mode, const std::vector<std::string>& lins, const std::string& out,
              const std::function<std::string(size_t)>& row) {
    std::vector<uint32_t> first_seen;
    std::vector<std::vector<uint32_t>> by(db.tids.size());
    std::vector<char> seen(db.tids.size(), 0);
    for (size_t i = 0; i < db.size(); ++i) {
        uint32_t t = db.e_tid[i];
        if (!seen[t]) { seen[t] = 1; first_seen.push_back(t); }
        by[t].push_back((uint32_t)i);
    }
    std::string ls;
    for (auto& l : lins) ls += (ls.empty() ? "" : ",") + l;
    write_text_file(out, [&](FILE* f) {
        fprintf(f, "#LQSNP\tmode=%s\tlineages=%s\tn=%zu\n", mode.c_str(), ls.c_str(), db.size());
        for (uint32_t t : first_seen) {
            auto& v = by[t];
            std::stable_sort(v.begin(), v.end(), [&](uint32_t a, uint32_t b) { return db.e_pos[a] < db.e_pos[b]; });
            for (uint32_t i : v) {
                std::string r = row(i);
                fprintf(f, "%s\t%u\t%s\n", db.tids[t].c_str(), db.e_pos[i], r.c_str());
            }
        }
    });
}
char one_base(char b) {
    b = (char)toupper((unsigned char)b);
    return (b == 'A' || b == 'C' || b == 'G' || b == 'T') ? b : '.';
}
}  // namespace

// ===================================================================================== snpdb-msa
namespace {
void msa_usage() {
    std::fprintf(stderr,
                 "usage: lq snpdb-msa --fasta COPIES.fa --cluster CLUSTERS --ancestry ANCESTRY --out SNP_DB_msa.lq.tsv\n"
                 "                              [--threads N] [--mafft-mode retree1|auto|fast] [--mafft PATH] [--tmpdir DIR]\n"
                 "                              [--legacy-lineages]\n"
                 "  per cluster: MAFFT alignment of the longest transcript of each copy; columns without gaps where the\n"
                 "  copies differ are SNPs; each copy's position gets the base expected for every lineage.\n"
                 "  Same database as build_snp_msa.py + export_snp_db.py of LineageQuant 1.0.0, with one fix:\n"
                 "  the lineage list is every lineage of the ancestry file. 1.0.0 took the lineages of the first SNP\n"
                 "  only (and looked the copies up by the id without its last \".\" field), which can leave lineages\n"
                 "  out of the database; --legacy-lineages reproduces that exactly.\n");
}
struct MsaCluster {
    std::string cid;
    std::vector<std::string> fids;   // full transcript ids (one per allele, insertion order)
};
struct MsaStats {
    long long clusters = 0, snps = 0, gaps = 0, identical = 0, failed = 0;
};
struct MsaEntry {
    std::string fid;
    uint32_t pos;
    std::string bases;   // per lineage of the universe, '\0' = none
};
}  // namespace

int snpdb_msa_main(int argc, char** argv) {
    std::string fasta, cluster, ancestry, out, mode = "retree1", mafft = "mafft", tmpbase;
    int threads = 0;
    bool legacy = false;
    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&]() -> std::string {
            if (i + 1 >= argc) { msa_usage(); die("missing value for " + a); }
            return argv[++i];
        };
        if (a == "--fasta") fasta = need();
        else if (a == "--cluster" || a == "--clusters") cluster = need();
        else if (a == "--ancestry") ancestry = need();
        else if (a == "--out" || a == "--lq-out") out = need();
        else if (a == "--threads" || a == "-t") threads = std::atoi(need().c_str());
        else if (a == "--mafft-mode" || a == "--mafft_mode") mode = need();
        else if (a == "--mafft") mafft = need();
        else if (a == "--tmpdir") tmpbase = need();
        else if (a == "--legacy-lineages") legacy = true;
        else if (a == "-h" || a == "--help") { msa_usage(); return 0; }
        else { msa_usage(); die("unknown argument " + a); }
    }
    if (fasta.empty() || cluster.empty() || ancestry.empty() || out.empty()) { msa_usage(); die("--fasta, --cluster, --ancestry and --out are required"); }
    if (threads <= 0) threads = std::max(1u, std::thread::hardware_concurrency() / 2);
    std::vector<std::string> mopts;
    if (mode == "auto") mopts = {"--auto", "--quiet"};
    else if (mode == "retree1") mopts = {"--retree", "1", "--quiet"};
    else if (mode == "fast") mopts = {"--retree", "1", "--maxiterate", "0", "--quiet"};
    else die("--mafft-mode must be auto, retree1 or fast");
    if (which(mafft).empty()) die("MAFFT not found on PATH (" + mafft + ")");

    // ancestry + lineage names
    std::unordered_map<std::string, std::string> anc;
    {
        LineReader rd(ancestry);
        std::string_view l;
        while (rd.next(l)) {
            if (py_strip(l).empty() || l.rfind("GeneID", 0) == 0) continue;
            auto p = py_split_ws(py_strip(l));
            if (p.size() >= 2) anc[std::string(before_dot(p[0]))] = upper(p[1]);
        }
    }
    std::set<std::string> linset;
    for (auto& kv : anc) linset.insert(kv.second);
    std::vector<std::string> lineages(linset.begin(), linset.end());
    // lineage universe: the ancestry lineages + UNKNOWN (copies whose lookup fails)
    std::vector<std::string> uni = lineages;
    if (!linset.count("UNKNOWN")) uni.push_back("UNKNOWN");
    std::unordered_map<std::string, int> uix;
    for (size_t i = 0; i < uni.size(); ++i) uix[uni[i]] = (int)i;
    std::string ls;
    for (auto& l : lineages) ls += (ls.empty() ? "" : ", ") + l;
    logm(sfmt("[INIT] MSA SNP | %d workers | %zu lineages: %s", threads, lineages.size(), ls.c_str()));

    // FASTA; best (longest, first on ties) transcript per gene
    std::unordered_map<std::string, std::string> seq;
    std::unordered_map<std::string, std::vector<std::string>> clean_to_full;
    std::vector<std::string> clean_order;
    read_fasta_bio(fasta, [&](std::string&& id, std::string&& s) {
        std::string c(before_dot(id));
        auto it = clean_to_full.find(c);
        if (it == clean_to_full.end()) { clean_order.push_back(c); clean_to_full[c].push_back(id); }
        else it->second.push_back(id);
        seq[id] = upper(s);
    });
    std::unordered_map<std::string, std::string> best;
    for (auto& c : clean_order) {
        auto& f = clean_to_full[c];
        size_t bi = 0;
        for (size_t i = 1; i < f.size(); ++i) if (seq[f[i]].size() > seq[f[bi]].size()) bi = i;
        best[c] = f[bi];
    }
    // clusters
    std::vector<MsaCluster> batches;
    long long dropped_seq = 0, dropped_anc = 0;
    {
        LineReader rd(cluster);
        std::string_view l;
        while (rd.next(l)) {
            if (py_strip(l).empty() || (!l.empty() && l[0] == '#')) continue;
            auto parts = py_split(py_strip(l), '\t');
            if (parts.size() < 2) die("cluster file: line without a second column: " + std::string(l.substr(0, 80)));
            MsaCluster C;
            C.cid = std::string(parts[0]);
            std::vector<std::string> keys;
            std::unordered_map<std::string, std::string> amap;
            for (auto a : py_split(parts[1], ',')) {
                std::string g(before_dot(a));
                auto b = best.find(g);
                if (b == best.end()) ++dropped_seq;
                else if (!anc.count(g)) ++dropped_anc;
                else if (amap.emplace(g, b->second).second) keys.push_back(g);
            }
            if (keys.size() >= 2) {
                for (auto& k : keys) C.fids.push_back(amap[k]);
                batches.push_back(std::move(C));
            }
        }
    }
    logm("[INFO] " + commas((long long)batches.size()) + " clusters to align");

    const size_t CH = 50;
    const size_t nchunks = (batches.size() + CH - 1) / CH;
    std::vector<std::vector<MsaEntry>> chunk_out(nchunks);
    std::vector<MsaStats> chunk_st(nchunks);
    std::atomic<size_t> next{0}, done{0};
    std::mutex logmx;
    const std::string troot = make_temp_dir(tmpbase, "snpdb_");
    auto worker = [&]() {
        const std::string wdir = make_temp_dir(troot, "mafft_");
        std::string cap;
        for (;;) {
            size_t c = next++;
            if (c >= nchunks) break;
            auto& E = chunk_out[c];
            auto& st = chunk_st[c];
            for (size_t bi = c * CH; bi < std::min(batches.size(), (c + 1) * CH); ++bi) {
                const MsaCluster& C = batches[bi];
                bool same = true;
                for (size_t i = 1; i < C.fids.size(); ++i) if (seq[C.fids[i]] != seq[C.fids[0]]) { same = false; break; }
                if (same) { ++st.identical; continue; }
                ++st.clusters;
                const std::string tf = wdir + "/c" + std::to_string(bi) + ".fa";
                {
                    FILE* f = fopen(tf.c_str(), "w");
                    if (!f) die("cannot write " + tf);
                    for (auto& id : C.fids) fprintf(f, ">%s\n%s\n", id.c_str(), seq[id].c_str());
                    if (fclose(f) != 0) die("error writing " + tf);
                }
                std::vector<std::string> cmd = {mafft};
                cmd.insert(cmd.end(), mopts.begin(), mopts.end());
                cmd.push_back(tf);
                cap.clear();
                RunResult r = run_proc(cmd, "|", "-", &cap, 120);
                unlink(tf.c_str());
                if (r.timed_out || r.status != 0 || py_strip(cap).empty()) { ++st.failed; continue; }
                // AlignIO.read(fasta): records of equal length
                std::vector<std::string> rid, rs;
                {
                    std::string_view all(cap);
                    bool in = false;
                    std::string id, s;
                    size_t b = 0;
                    auto flush = [&]() {
                        std::string t;
                        for (char ch : s) if (ch != ' ' && ch != '\r') t += ch;
                        rid.push_back(id);
                        rs.push_back(t);
                        s.clear();
                    };
                    while (b <= all.size()) {
                        size_t e = all.find('\n', b);
                        std::string_view ln = all.substr(b, e == std::string_view::npos ? std::string_view::npos : e - b);
                        if (!ln.empty() && ln.back() == '\r') ln.remove_suffix(1);
                        if (!ln.empty() && ln[0] == '>') {
                            if (in) flush();
                            in = true;
                            auto w = py_split_ws(py_rstrip(ln.substr(1)));
                            id = w.empty() ? std::string() : std::string(w[0]);
                        } else if (in) { auto rr = py_rstrip(ln); s.append(rr.data(), rr.size()); }
                        if (e == std::string_view::npos) break;
                        b = e + 1;
                    }
                    if (in) flush();
                }
                bool ok = !rs.empty();
                for (auto& x : rs) if (x.size() != rs[0].size()) ok = false;
                if (!ok) { ++st.failed; continue; }
                const size_t R = rs.size(), L = rs[0].size();
                std::vector<int> lin(R);
                for (size_t i = 0; i < R; ++i) {
                    std::string k = legacy ? std::string(rsplit_dot(rid[i])) : std::string(before_dot(rid[i]));
                    auto a = anc.find(k);
                    lin[i] = uix[a == anc.end() ? std::string("UNKNOWN") : a->second];
                }
                std::vector<uint32_t> up(R, 0);
                for (size_t ci = 0; ci < L; ++ci) {
                    bool gap = false, diff = false;
                    for (size_t i = 0; i < R; ++i) {
                        char ch = rs[i][ci];
                        if (ch == '-' || ch == 'N') gap = true;
                        if (ch != rs[0][ci]) diff = true;
                    }
                    if (gap) {
                        ++st.gaps;
                        for (size_t i = 0; i < R; ++i) if (rs[i][ci] != '-') ++up[i];
                        continue;
                    }
                    if (diff) {
                        std::string bases(uni.size(), '\0');
                        for (size_t i = 0; i < R; ++i) if (!bases[lin[i]]) bases[lin[i]] = rs[i][ci];
                        for (size_t i = 0; i < R; ++i) E.push_back(MsaEntry{rid[i], up[i], bases});
                        ++st.snps;
                    }
                    for (size_t i = 0; i < R; ++i) ++up[i];   // no '-' in this column
                }
            }
            size_t d = ++done;
            if (d % 50 == 0) {
                std::lock_guard<std::mutex> g(logmx);
                logm("  ... " + std::to_string(d) + "/" + std::to_string(nchunks) + " chunks");
            }
        }
        remove_dir_tree(wdir);
    };
    std::vector<std::thread> ws;
    for (int t = 0; t < threads; ++t) ws.emplace_back(worker);
    for (auto& t : ws) t.join();
    remove_dir_tree(troot);

    // merge in chunk order (= 1.0.0 with one worker; with several workers 1.0.0 merged in completion order)
    SnpDb db;
    MsaStats agg;
    for (size_t c = 0; c < nchunks; ++c) {
        agg.clusters += chunk_st[c].clusters; agg.snps += chunk_st[c].snps; agg.gaps += chunk_st[c].gaps;
        agg.identical += chunk_st[c].identical; agg.failed += chunk_st[c].failed;
        for (auto& e : chunk_out[c]) db.put(db.tid_of(e.fid), e.pos, std::move(e.bases));
        std::vector<MsaEntry>().swap(chunk_out[c]);
    }
    logm(sfmt("[SUMMARY] clusters=%s identical=%s SNPs=%s DB_entries=%s%s", commas(agg.clusters).c_str(),
              commas(agg.identical).c_str(), commas(agg.snps).c_str(), commas((long long)db.size()).c_str(),
              agg.failed ? sfmt(" (MAFFT failed / timed out on %s clusters)", commas(agg.failed).c_str()).c_str() : ""));
    std::vector<std::string> lins;
    if (legacy) {
        if (db.size() == 0) die("could not infer lineage names from the SNP DB");
        const std::string& b = db.e_val[0];
        for (size_t i = 0; i < uni.size(); ++i) if (b[i]) lins.push_back(uni[i]);
        std::sort(lins.begin(), lins.end());
    } else lins = lineages;
    std::vector<int> col;
    for (auto& l : lins) col.push_back(uix.count(l) ? uix[l] : -1);
    write_lq(db, "msa", lins, out, [&](size_t i) {
        std::string r;
        for (int c : col) r += (c >= 0 && db.e_val[i][c]) ? one_base(db.e_val[i][c]) : '.';
        return r;
    });
    std::string lj;
    for (auto& l : lins) lj += (lj.empty() ? "" : ", ") + l;
    logm("[SUCCESS] " + commas((long long)db.size()) + " SNPs on " + commas((long long)std::set<uint32_t>(db.e_tid.begin(), db.e_tid.end()).size()) +
         " transcripts, lineages [" + lj + "] -> " + out);
    if (dropped_seq || dropped_anc)
        logm("[INFO] cluster members skipped: " + commas(dropped_seq) + " without a sequence, " + commas(dropped_anc) + " without a lineage");
    return 0;
}

// ===================================================================================== snpdb-vcf
namespace {
void vcf_usage() {
    std::fprintf(stderr,
                 "usage: lq snpdb-vcf --vcf POP.vcf[.gz] --gff3 REF.gff3 --sample-lists NAME:FILE [NAME:FILE ...]\n"
                 "                              --out SNP_DB_population.lq.tsv [--min_dp 3] [--min_called 3] [--max_dp 200]\n"
                 "                              [--min_af_diff 0]\n"
                 "  biallelic SNVs projected from genome to transcript coordinates (GFF3 CDS/exon features);\n"
                 "  per lineage: allele frequency = alt / (ref + alt) of the AD of its samples with min_dp <= DP <= max_dp\n"
                 "  (0.5 with fewer than min_called such samples). Same database as build_snp_vcf.py + export_snp_db.py\n"
                 "  --mode population of LineageQuant 1.0.0. Plain or gzip/bgzip VCF (not BCF).\n");
}
struct Ival {
    int64_t s, e, off;
    uint32_t key, rank, tx;
    bool minus;
};
const int32_t I32_MISSING = INT32_MIN, I32_VEND = INT32_MIN + 1;
// one integer FORMAT value as htslib stores it
int32_t fmt_int(std::string_view v, bool& ok) {
    ok = true;
    if (v.empty() || v == ".") return I32_MISSING;
    long long x;
    if (!py_int(v, x)) { ok = false; return I32_MISSING; }
    return (int32_t)x;
}
}  // namespace

int snpdb_vcf_main(int argc, char** argv) {
    std::string vcf, gff, out;
    std::vector<std::string> lists;
    long long min_dp = 3, min_called = 3, max_dp = 200;
    double min_af_diff = 0.0;
    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&]() -> std::string {
            if (i + 1 >= argc) { vcf_usage(); die("missing value for " + a); }
            return argv[++i];
        };
        if (a == "--vcf") vcf = need();
        else if (a == "--gff3" || a == "--gff") gff = need();
        else if (a == "--out" || a == "--lq-out") out = need();
        else if (a == "--sample-lists") {
            while (i + 1 < argc && argv[i + 1][0] != '-') lists.push_back(argv[++i]);
        }
        else if (a == "--min_dp" || a == "--min-dp") min_dp = std::atoll(need().c_str());
        else if (a == "--min_called" || a == "--min-called") min_called = std::atoll(need().c_str());
        else if (a == "--max_dp" || a == "--max-dp") max_dp = std::atoll(need().c_str());
        else if (a == "--min_af_diff" || a == "--min-af-diff") { if (!py_float(need(), min_af_diff)) die("--min_af_diff: not a number"); }
        else if (a == "-h" || a == "--help") { vcf_usage(); return 0; }
        else { vcf_usage(); die("unknown argument " + a); }
    }
    if (vcf.empty() || gff.empty() || out.empty() || lists.empty()) { vcf_usage(); die("--vcf, --gff3, --sample-lists and --out are required"); }
    const double T0 = wall();
    std::map<std::string, std::unordered_set<std::string>> groups;   // sorted names
    for (auto& it : lists) {
        size_t c = it.find(':');
        if (c == std::string::npos) die("--sample-lists must be NAME:FILE, got " + it);
        std::unordered_set<std::string> s;
        LineReader rd(it.substr(c + 1));
        std::string_view l;
        while (rd.next(l)) { auto x = py_strip(l); if (!x.empty()) s.insert(std::string(x)); }
        groups[it.substr(0, c)] = std::move(s);
    }
    std::vector<std::string> lins;
    for (auto& g : groups) lins.push_back(g.first);
    {
        std::string s;
        for (auto& l : lins) s += (s.empty() ? "" : ", ") + l;
        logm("[INFO] " + std::to_string(lins.size()) + " lineages: [" + s + "]");
    }

    // ---- genome -> transcript coordinates (build_genomic_to_transcript_map, as intervals)
    logm("[STEP 1] Projecting genomic->transcript coordinates...");
    struct TxKey { std::string chrom, tx, strand; };
    std::vector<TxKey> keys;
    std::map<std::tuple<std::string, std::string, std::string>, uint32_t> kix;
    std::vector<std::vector<std::pair<int64_t, int64_t>>> exons;
    {
        LineReader rd(gff);
        std::string_view l;
        while (rd.next(l)) {
            if ((!l.empty() && l[0] == '#') || py_strip(l).empty()) continue;
            auto p = py_split(py_strip(l), '\t');
            if (p.size() < 9) continue;
            std::string ty = lower_ascii(p[2]);
            if (ty != "cds" && ty != "exon") continue;
            long long st, en;
            if (!py_int(p[3], st) || !py_int(p[4], en)) die("GFF: bad coordinates: " + std::string(l.substr(0, 100)));
            std::string tx;
            bool found = false;
            for (auto kv : py_split(p[8], ';')) {
                auto k = py_strip(kv);
                if (k.rfind("Parent=", 0) == 0) {
                    auto f = py_split(k, '=');
                    tx = std::string(before_dot(f[1]));
                    found = true;
                    break;
                }
            }
            if (!found || tx.empty()) continue;
            auto key = std::make_tuple(std::string(p[0]), tx, std::string(p[6]));
            auto it = kix.find(key);
            uint32_t id;
            if (it == kix.end()) {
                id = (uint32_t)keys.size();
                kix.emplace(key, id);
                keys.push_back({std::string(p[0]), tx, std::string(p[6])});
                exons.emplace_back();
            } else id = it->second;
            exons[id].emplace_back(st, en);
        }
    }
    std::unordered_map<std::string, std::vector<Ival>> ivs;   // chrom -> intervals
    std::vector<std::string> txnames;
    std::unordered_map<std::string, uint32_t> txix;
    long long npos = 0;
    for (uint32_t k = 0; k < keys.size(); ++k) {
        auto ex = exons[k];
        const bool minus = keys[k].strand == "-";
        if (minus) std::stable_sort(ex.begin(), ex.end(), [](auto& a, auto& b) { return a.first > b.first; });
        else std::stable_sort(ex.begin(), ex.end(), [](auto& a, auto& b) { return a.first < b.first; });
        auto t = txix.find(keys[k].tx);
        uint32_t tx;
        if (t == txix.end()) { tx = (uint32_t)txnames.size(); txix.emplace(keys[k].tx, tx); txnames.push_back(keys[k].tx); }
        else tx = t->second;
        int64_t off = 0;
        uint32_t rank = 0;
        for (auto& e : ex) {
            if (e.second >= e.first) {
                ivs[keys[k].chrom].push_back(Ival{e.first, e.second, off, k, rank, tx, minus});
                npos += e.second - e.first + 1;
            }
            ++rank;
            off += e.second - e.first + 1;
        }
    }
    // bins of 64 kb -> intervals overlapping the bin
    const int SH = 16;
    std::unordered_map<std::string, std::unordered_map<int64_t, std::vector<uint32_t>>> bins;
    for (auto& kv : ivs) {
        auto& B = bins[kv.first];
        for (uint32_t i = 0; i < kv.second.size(); ++i)
            for (int64_t b = kv.second[i].s >> SH; b <= (kv.second[i].e >> SH); ++b) B[b].push_back(i);
    }
    logm("  " + commas(npos) + " genomic CDS/exon base positions (with overlaps)");

    // ---- VCF
    logm("[STEP 2] Building SNP DB from VCF...");
    LineReader rd(vcf);
    std::string_view l;
    std::vector<std::vector<int>> idx(lins.size());
    bool header = false;
    long long n_total = 0, n_cds = 0, n_entries = 0;
    SnpDb db;
    std::vector<int32_t> dp, ad;
    std::vector<std::pair<std::pair<uint32_t, uint32_t>, std::pair<uint32_t, int64_t>>> hits;   // ((key, rank), (tx, tpos))
    while (rd.next(l)) {
        if (!header) {
            if (l.rfind("BCF", 0) == 0) die(vcf + " looks like BCF; convert it with: bcftools view -Oz -o x.vcf.gz " + vcf);
            if (l.rfind("##", 0) == 0) continue;
            if (l.rfind("#CHROM", 0) == 0) {
                auto h = split_tab(l);
                std::vector<std::string> samples;
                for (size_t i = 9; i < h.size(); ++i) samples.emplace_back(h[i]);
                for (size_t li = 0; li < lins.size(); ++li) {
                    for (size_t s = 0; s < samples.size(); ++s) if (groups[lins[li]].count(samples[s])) idx[li].push_back((int)s);
                    if (idx[li].empty()) logm("[WARNING] lineage '" + lins[li] + "' has no matching samples in VCF header.");
                }
                header = true;
                continue;
            }
            die(vcf + ": no #CHROM header line before the records");
        }
        if (l.empty()) continue;
        auto f = split_tab(l);
        if (f.size() < 8) die(vcf + ": record with fewer than 8 columns");
        std::string_view ref = f[3], alt = f[4];
        if (alt == "." || alt.find(',') != std::string_view::npos || ref.size() > 1 || alt.size() > 1) continue;
        ++n_total;
        long long pos;
        if (!py_int(f[1], pos)) die(vcf + ": bad POS " + std::string(f[1]));
        auto bc = bins.find(std::string(f[0]));
        if (bc == bins.end()) continue;
        auto bb = bc->second.find(pos >> SH);
        if (bb == bc->second.end()) continue;
        hits.clear();
        const auto& IV = ivs[std::string(f[0])];
        for (uint32_t i : bb->second) {
            const Ival& v = IV[i];
            if (pos < v.s || pos > v.e) continue;
            int64_t tpos = v.off + (v.minus ? v.e - pos : pos - v.s);
            hits.push_back({{v.key, v.rank}, {v.tx, tpos}});
        }
        if (hits.empty()) continue;
        std::sort(hits.begin(), hits.end(), [](auto& a, auto& b) { return a.first < b.first; });
        // FORMAT: AD and DP (cyvcf2 rec.format() is None when the record lacks the field)
        if (f.size() < 9) continue;
        int iad = -1, idp = -1;
        {
            auto fm = py_split(f[8], ':');
            for (size_t i = 0; i < fm.size(); ++i) { if (fm[i] == "AD") iad = (int)i; else if (fm[i] == "DP") idp = (int)i; }
        }
        if (iad < 0 || idp < 0) continue;
        const size_t NS = f.size() - 9;
        dp.assign(NS, I32_MISSING);
        std::vector<std::vector<int32_t>> adv(NS);
        size_t nmax = 0;
        for (size_t s = 0; s < NS; ++s) {
            auto sf = py_split(f[9 + s], ':');
            bool ok;
            if ((int)sf.size() > idp) { dp[s] = fmt_int(sf[idp], ok); if (!ok) die(vcf + ": bad DP value at " + std::string(f[0]) + ":" + std::string(f[1])); }
            if ((int)sf.size() > iad && sf[iad] != "." && !sf[iad].empty()) {
                for (auto x : py_split(sf[iad], ',')) { adv[s].push_back(fmt_int(x, ok)); if (!ok) die(vcf + ": bad AD value at " + std::string(f[0]) + ":" + std::string(f[1])); }
            } else adv[s].push_back(I32_MISSING);
            nmax = std::max(nmax, adv[s].size());
        }
        std::vector<std::string> afs(lins.size());
        std::vector<double> afv(lins.size());
        for (size_t li = 0; li < lins.size(); ++li) {
            double af = 0.5;
            if (!idx[li].empty()) {
                long long n_ok = 0, sa = 0, sr = 0;
                for (int s : idx[li]) if (dp[s] >= min_dp && dp[s] <= max_dp) ++n_ok;
                if (n_ok >= min_called) {
                    if (nmax < 2) die(vcf + ": AD with one value per sample at " + std::string(f[0]) + ":" + std::string(f[1]) + " (needs REF,ALT)");
                    for (int s : idx[li]) {
                        if (!(dp[s] >= min_dp && dp[s] <= max_dp)) continue;
                        auto& v = adv[s];
                        sr += v.size() > 0 ? v[0] : I32_VEND;
                        sa += v.size() > 1 ? v[1] : I32_VEND;
                    }
                    long long tot = sa + sr;
                    af = tot > 0 ? (double)sa / (double)tot : 0.5;
                }
            }
            afv[li] = af;
        }
        if (min_af_diff > 0 && !afv.empty()) {
            double mx = *std::max_element(afv.begin(), afv.end()), mn = *std::min_element(afv.begin(), afv.end());
            if (mx - mn < min_af_diff) continue;
        }
        ++n_cds;
        std::string val;
        val += ref.size() == 1 ? (char)toupper((unsigned char)ref[0]) : '.';
        val += '\t';
        val += alt.size() == 1 ? (char)toupper((unsigned char)alt[0]) : '.';
        val += '\t';
        for (size_t li = 0; li < lins.size(); ++li) { if (li) val += ','; val += py_repr(afv[li]); }
        for (auto& h : hits) {
            if (h.second.second < 0 || h.second.second > (int64_t)UINT32_MAX) continue;
            db.put(db.tid_of(txnames[h.second.first]), (uint32_t)h.second.second, std::string(val));
            ++n_entries;
        }
    }
    logm(sfmt("[SUMMARY] biallelic=%s in_CDS=%s entries=%s", commas(n_total).c_str(), commas(n_cds).c_str(), commas(n_entries).c_str()));
    if (db.size() == 0 || lins.empty()) die("could not infer lineage names from the SNP DB (no SNP in the transcripts)");
    write_lq(db, "population", lins, out, [&](size_t i) { return db.e_val[i]; });
    logm("[SUCCESS] " + commas((long long)db.size()) + " SNPs on " + commas((long long)std::set<uint32_t>(db.e_tid.begin(), db.e_tid.end()).size()) +
         " transcripts -> " + out + sfmt(" (%.0fs)", wall() - T0));
    return 0;
}

// ===================================================================================== build-markers
namespace {
void markers_usage() {
    std::fprintf(stderr,
                 "usage: lq build-markers --lineage NAME:FASTA [--lineage ...] -o OUTDIR [-k 31] [-t 8]\n"
                 "                                  [--polyploid NAME:FASTA ...] [--mode ancestor|ancestor_x_poly|poly]\n"
                 "                                  [--at-gc-frac 0.8] [--keep-intermediate]\n"
                 "  marker k-mers of lineage L (KMC, canonical k-mers):\n"
                 "    ancestor         (default without --polyploid; what `lq run` uses)  anc_L - other anc\n"
                 "    ancestor_x_poly  (default with --polyploid)  (anc_L - other anc) * poly_L - other poly\n"
                 "    poly             poly_L - other poly\n"
                 "  then low-complexity k-mers are dropped; OUTDIR/<NAME>.fasta per lineage. Same as build_markers.py.\n"
                 "  The same NAME given twice = several inputs of one lineage (concatenated).\n");
}
}  // namespace

bool is_low_complexity(const std::string& km, int k, int thr) {
    for (const char* h : {"AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG"}) if (km.find(h) != std::string::npos) return true;
    int a = 0, t = 0, c = 0, g = 0;
    for (char x : km) { a += x == 'A'; t += x == 'T'; c += x == 'C'; g += x == 'G'; }
    if (a + t > thr || c + g > thr) return true;
    bool one = !km.empty();
    for (char x : km) if (x != km[0]) { one = false; break; }
    if (one) return true;
    {
        std::string di = km.substr(0, 2), rep;
        for (int i = 0; i < k / 2; ++i) rep += di;
        if (rep == km.substr(0, std::min(km.size(), (size_t)((k / 2) * 2)))) return true;
    }
    // re.search(r'([ACGT]{2,10})\1{2,}', kmer)
    const size_t n = km.size();
    auto acgt = [](char x) { return x == 'A' || x == 'C' || x == 'G' || x == 'T'; };
    for (size_t i = 0; i < n; ++i) {
        for (size_t L = 2; L <= 10 && i + 3 * L <= n; ++L) {
            bool ok = true;
            for (size_t j = 0; j < L; ++j) if (!acgt(km[i + j])) { ok = false; break; }
            if (!ok) break;   // a longer motif contains the same non-ACGT base
            if (km.compare(i + L, L, km, i, L) == 0 && km.compare(i + 2 * L, L, km, i, L) == 0) return true;
        }
    }
    return false;
}

std::map<std::string, std::string> build_markers(const std::vector<std::pair<std::string, std::string>>& anc,
                                                 const std::string& outdir, int k, int threads, double at_gc_frac,
                                                 const std::vector<std::pair<std::string, std::string>>& poly,
                                                 const std::string& mode, bool keep_intermediate,
                                                 std::vector<std::string>* order) {
    for (const char* t : {"kmc", "kmc_tools"}) if (which(t).empty()) die(std::string("Required tool '") + t + "' not found on PATH.");
    mkdirs(outdir);
    std::vector<std::string> names;
    std::map<std::string, std::string> A, P;
    for (auto& x : anc) A[x.first] = x.second;
    for (auto& x : poly) P[x.first] = x.second;
    if (mode == "poly") {
        if (poly.size() < 2) die("mode=poly needs >= 2 --polyploid");
        for (auto& x : poly) names.push_back(x.first);
    } else {
        if (anc.size() < 2) die("need >= 2 lineages");
        for (auto& x : anc) names.push_back(x.first);
        if (mode == "ancestor_x_poly")
            for (auto& n : names) if (!P.count(n)) die("mode=ancestor_x_poly: no --polyploid for " + n);
    }
    auto kmc_db = [&](const std::string& fa, const std::string& prefix) {
        if (file_exists(prefix + ".kmc_pre")) return prefix;
        const std::string tmp = prefix + "_kmctmp";
        mkdirs(tmp);
        RunResult r = run_proc({"kmc", "-k" + std::to_string(k), "-ci1", "-cs1000000", "-t" + std::to_string(threads), "-fm", fa, prefix, tmp}, "-", "-");
        remove_dir_tree(tmp);
        if (r.status != 0) die("kmc failed (exit " + std::to_string(r.status) + ") on " + fa);
        return prefix;
    };
    std::map<std::string, std::string> adb, pdb;
    if (mode != "poly")
        for (auto& n : names) { logm("[build_markers] KMC count ancestor: " + n); adb[n] = kmc_db(A[n], outdir + "/anc_" + n + "_k" + std::to_string(k)); }
    if (mode != "ancestor")
        for (auto& n : names) { logm("[build_markers] KMC count polyploid subgenome: " + n); pdb[n] = kmc_db(P[n], outdir + "/poly_" + n + "_k" + std::to_string(k)); }

    std::map<std::string, std::string> txts;
    for (auto& n : names) {
        std::vector<std::string> others;
        for (auto& o : names) if (o != n) others.push_back(o);
        const std::string spec = outdir + "/" + n + "_spec", txt = outdir + "/" + n + ".txt";
        std::vector<std::string> lines = {"INPUT:"};
        std::string expr;
        if (mode == "poly") {
            lines.push_back("s = " + pdb[n] + " -ci1");
            for (size_t i = 0; i < others.size(); ++i) lines.push_back("q" + std::to_string(i + 1) + " = " + pdb[others[i]] + " -ci1");
            expr = "s";
            for (size_t i = 0; i < others.size(); ++i) expr += " - q" + std::to_string(i + 1);
        } else {
            lines.push_back("s = " + adb[n] + " -ci1");
            for (size_t i = 0; i < others.size(); ++i) lines.push_back("o" + std::to_string(i + 1) + " = " + adb[others[i]] + " -ci1");
            expr = "s";
            for (size_t i = 0; i < others.size(); ++i) expr += " - o" + std::to_string(i + 1);
            if (mode == "ancestor_x_poly") {
                lines.push_back("p = " + pdb[n] + " -ci1");
                for (size_t i = 0; i < others.size(); ++i) lines.push_back("q" + std::to_string(i + 1) + " = " + pdb[others[i]] + " -ci1");
                std::string e2 = "(" + expr + ") * p";
                for (size_t i = 0; i < others.size(); ++i) e2 += " - q" + std::to_string(i + 1);
                expr = e2;
            }
        }
        lines.push_back("OUTPUT:");
        lines.push_back(spec + " = " + expr);
        logm("[build_markers] " + n + ": " + expr);
        const std::string cfg = outdir + "/" + n + ".cfg";
        write_text_file(cfg, [&](FILE* f) {
            for (size_t i = 0; i < lines.size(); ++i) { fputs(lines[i].c_str(), f); fputc('\n', f); }
        });
        run_checked({"kmc_tools", "-t" + std::to_string(threads), "complex", cfg}, "-", "-");
        run_checked({"kmc_tools", "transform", spec, "dump", txt}, "-", "-");
        txts[n] = txt;
    }
    // merged table -> low-complexity filter -> one FASTA per lineage (in order of first kept k-mer)
    const int thr = (int)(at_gc_frac * (double)k);
    std::map<std::string, std::string> paths;
    std::map<std::string, FILE*> fh;
    std::map<std::string, std::pair<long long, long long>> kr;   // kept, removed
    FILE* fm = keep_intermediate ? fopen((outdir + "/merged.txt").c_str(), "w") : nullptr;
    FILE* fl = keep_intermediate ? fopen((outdir + "/merged_lc.txt").c_str(), "w") : nullptr;
    if (fm) fputs("Kmer\tFinalScore\tSpecies\n", fm);
    if (fl) fputs("Kmer\tFinalScore\tSpecies\n", fl);
    for (auto& n : names) {
        LineReader rd(txts[n]);
        std::string_view l;
        while (rd.next(l)) {
            auto f0 = l.substr(0, l.find('\t'));
            auto w = py_split_ws(f0);
            if (w.empty()) die(txts[n] + ": empty k-mer line");
            std::string km(w[0]);
            if (km.empty()) continue;
            if (fm) fprintf(fm, "%s\t1\t%s\n", km.c_str(), n.c_str());
            if (is_low_complexity(km, k, thr)) { kr[n].second++; continue; }
            kr[n].first++;
            if (fl) fprintf(fl, "%s\t1\t%s\n", km.c_str(), n.c_str());
            auto it = fh.find(n);
            if (it == fh.end()) {
                std::string p = outdir + "/" + n + ".fasta";
                FILE* f = fopen(p.c_str(), "w");
                if (!f) die("cannot write " + p);
                it = fh.emplace(n, f).first;
                paths[n] = p;
                if (order) order->push_back(n);
            }
            fprintf(it->second, ">%s-1\n%s\n", n.c_str(), km.c_str());
        }
    }
    for (auto& kv : fh) if (fclose(kv.second) != 0) die("error writing " + paths[kv.first]);
    if (fm) fclose(fm);
    if (fl) fclose(fl);
    long long tk = 0, tr = 0;
    logm("[filter_low_complexity] per-lineage kept / removed:");
    for (auto& kv : kr) {
        long long a = kv.second.first, r = kv.second.second;
        tk += a; tr += r;
        logm(sfmt("  %-16s kept=%10lld  removed=%10lld  (%.1f%% LC)", kv.first.c_str(), a, r, a + r ? 100.0 * r / (a + r) : 0.0));
    }
    logm(sfmt("  %-16s kept=%10lld  removed=%10lld  (%.1f%%)", "TOTAL", tk, tr, tk + tr ? 100.0 * tr / (tk + tr) : 0.0));
    if (!keep_intermediate) for (auto& kv : txts) unlink(kv.second.c_str());
    logm("[build_markers] final marker counts:");
    for (auto& n : names) logm("    " + n + ": " + commas(kr[n].first));
    return paths;
}

int build_markers_main(int argc, char** argv) {
    std::vector<std::pair<std::string, std::vector<std::string>>> lin, pol;
    std::string outdir, mode;
    int k = 31, threads = 8;
    double frac = 0.8;
    bool keep = false;
    auto add = [](std::vector<std::pair<std::string, std::vector<std::string>>>& g, const std::string& it, const char* flag) {
        size_t c = it.find(':');
        if (c == std::string::npos) die(std::string(flag) + " must be NAME:PATH, got: " + it);
        std::string n = it.substr(0, c);
        for (auto& x : g) if (x.first == n) { x.second.push_back(it.substr(c + 1)); return; }
        g.push_back({n, {it.substr(c + 1)}});
    };
    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&]() -> std::string {
            if (i + 1 >= argc) { markers_usage(); die("missing value for " + a); }
            return argv[++i];
        };
        if (a == "--lineage") add(lin, need(), "--lineage");
        else if (a == "--polyploid") add(pol, need(), "--polyploid");
        else if (a == "--mode") mode = need();
        else if (a == "-o" || a == "--outdir") outdir = need();
        else if (a == "-k" || a == "--kmer") k = std::atoi(need().c_str());
        else if (a == "-t" || a == "--threads") threads = std::atoi(need().c_str());
        else if (a == "--at-gc-frac") { if (!py_float(need(), frac)) die("--at-gc-frac: not a number"); }
        else if (a == "--keep-intermediate") keep = true;
        else if (a == "-h" || a == "--help") { markers_usage(); return 0; }
        else { markers_usage(); die("unknown argument " + a); }
    }
    if (outdir.empty()) { markers_usage(); die("-o OUTDIR is required"); }
    if (mode.empty()) mode = pol.empty() ? "ancestor" : "ancestor_x_poly";
    if (mode != "ancestor" && mode != "ancestor_x_poly" && mode != "poly") die("--mode must be ancestor, ancestor_x_poly or poly");
    auto merge = [&](const std::vector<std::pair<std::string, std::vector<std::string>>>& g, const std::string& tag) {
        std::vector<std::pair<std::string, std::string>> r;
        for (auto& x : g) {
            if (x.second.size() == 1) { r.push_back({x.first, x.second[0]}); continue; }
            mkdirs(outdir);
            std::string m = outdir + "/_" + tag + "_" + x.first + "_merged.fa";
            concat_files(x.second, m);
            logm("[build_markers] " + tag + " " + x.first + ": merged " + std::to_string(x.second.size()) + " inputs");
            r.push_back({x.first, m});
        }
        return r;
    };
    auto A = merge(lin, "anc"), P = merge(pol, "poly");
    if (mode != "poly" && A.size() < 2) die("--lineage: at least two lineages");
    if (mode != "ancestor" && P.size() < 2) die("mode=" + mode + " needs at least two --polyploid");
    std::vector<std::string> order;
    auto paths = build_markers(A, outdir, k, threads, frac, P, mode, keep, &order);
    std::printf("[done] marker FASTAs:\n");
    for (auto& n : order) std::printf("  %s: %s\n", n.c_str(), paths[n].c_str());
    return 0;
}

void concat_files(const std::vector<std::string>& in, const std::string& out) {
    FILE* o = fopen(out.c_str(), "wb");
    if (!o) die("cannot write " + out);
    std::vector<char> buf(1 << 20);
    for (auto& p : in) {
        FILE* f = fopen(p.c_str(), "rb");
        if (!f) die("cannot open " + p);
        size_t n;
        while ((n = fread(buf.data(), 1, buf.size(), f)) > 0)
            if (fwrite(buf.data(), 1, n, o) != n) die("error writing " + out);
        fclose(f);
    }
    if (fclose(o) != 0) die("error writing " + out);
}

// ===================================================================================== extract-transcripts (genome_pair.py)
namespace {
// Python text-mode line iteration (universal newlines): lines with their "\n"; "\r\n" and "\r" become "\n"
void read_text_lines(const std::string& path, const std::function<void(const std::string&)>& fn) {
    gzFile g = gzopen(path.c_str(), "rb");
    if (!g) die("cannot open " + path);
    gzbuffer(g, 1 << 20);
    std::vector<char> buf(1 << 20);
    std::string cur;
    bool cr = false;
    int n;
    while ((n = gzread(g, buf.data(), (unsigned)buf.size())) > 0) {
        for (int i = 0; i < n; ++i) {
            char c = buf[i];
            if (cr) { cr = false; if (c == '\n') continue; }
            if (c == '\r') { cur += '\n'; fn(cur); cur.clear(); cr = true; continue; }
            cur += c;
            if (c == '\n') { fn(cur); cur.clear(); }
        }
    }
    gzclose(g);
    if (!cur.empty()) fn(cur);
}
std::string norm_token(const std::string& s0) {
    std::string s;
    for (char c : lower_ascii(s0)) if ((c >= 'a' && c <= 'z') || (c >= '0' && c <= '9')) s += c;
    for (const char* p : {"chromosome", "chrom", "chr", "scaffold", "contig", "seq"})
        if (s.rfind(p, 0) == 0) return s.substr(strlen(p));
    return s;
}
}  // namespace

// samtools faidx (htslib fai_build_core): NAME LENGTH OFFSET LINEBASES LINEWIDTH; false + message on a format error
bool write_fai(const std::string& fa, const std::string& fai, std::vector<std::pair<std::string, long long>>* lens, std::string* err) {
    FILE* f = fopen(fa.c_str(), "rb");
    if (!f) { *err = "cannot open " + fa; return false; }
    std::vector<char> buf(1 << 20);
    size_t bn = 0, bp = 0;
    uint64_t off = 0;
    auto getc_ = [&]() -> int {
        if (bp == bn) { bn = fread(buf.data(), 1, buf.size(), f); bp = 0; if (bn == 0) return -1; }
        ++off;
        return (unsigned char)buf[bp++];
    };
    struct Rec { std::string name; long long len; uint64_t off; int lb, lw; };
    std::vector<Rec> recs;
    std::unordered_set<std::string> names;
    enum { OUT, NAME, SEQ } state = OUT;
    bool read_done = false;
    std::string name;
    long long seq_len = 0, line_len = 0, char_len = 0;
    uint64_t seq_off = 0;
    int line_num = 1, c;
    auto insert = [&]() {
        if (!names.insert(name).second) {
            fprintf(stderr, "[W::fai_build_core] Ignoring duplicate sequence \"%s\" at byte offset %llu\n", name.c_str(), (unsigned long long)seq_off);
            return;
        }
        recs.push_back({name, seq_len, seq_off, (int)char_len, (int)line_len});
    };
    bool fail = false;
    while (!fail && (c = getc_()) >= 0) {
        switch (state) {
        case OUT:
            if (c == '>') state = NAME;
            else if (c == '\r') { if ((c = getc_()) == '\n') ++line_num; else { *err = "Format error, carriage return not followed by new line at line " + std::to_string(line_num); fail = true; } }
            else if (c == '\n') ++line_num;
            else { *err = "Format error, unexpected character at line " + std::to_string(line_num); fail = true; }
            break;
        case NAME:
            if (read_done) { insert(); read_done = false; }
            name.clear();
            do {
                if (!isspace(c)) name += (char)c;
                else if (!name.empty() || c == '\n') break;
            } while ((c = getc_()) >= 0);
            if (c < 0) { *err = "The last entry '" + name + "' has no sequence"; fail = true; break; }
            if (c != '\n') while ((c = getc_()) >= 0 && c != '\n') {}
            state = SEQ;
            seq_len = line_len = char_len = 0;
            seq_off = off;
            ++line_num;
            break;
        case SEQ: {
            if (c == '\n') { state = OUT; ++line_num; continue; }
            if (c == '>') { state = NAME; continue; }
            long long ll = 0, cl = 0;
            read_done = true;
            do {
                ++ll;
                if (isgraph(c)) ++cl;
            } while ((c = getc_()) >= 0 && c != '\n');
            ++ll;
            seq_len += cl;
            if (line_len == 0) { line_len = ll; char_len = cl; }
            else if (line_len > ll) state = OUT;
            else if (line_len < ll) { *err = "Different line length in sequence '" + name + "'"; fail = true; break; }
            ++line_num;
            break;
        }
        }
    }
    fclose(f);
    if (fail) return false;
    if (read_done) insert();
    write_text_file(fai, [&](FILE* o) {
        for (auto& r : recs) fprintf(o, "%s\t%lld\t%llu\t%d\t%d\n", r.name.c_str(), r.len, (unsigned long long)r.off, r.lb, r.lw);
    });
    if (lens) for (auto& r : recs) lens->push_back({r.name, r.len});
    return true;
}

void genome_pair_extract(const std::string& fasta, const std::string& gff, const std::string& workdir0, const std::string& name,
                         const std::string& out_fa, const std::string& id_prefix, long long min_chr_len) {
    if (!file_exists(fasta)) die("FASTA not found: " + fasta);
    if (!file_exists(gff)) die("GFF not found: " + gff);
    if (which("gffread").empty()) die("Required tool 'gffread' not found on PATH.");
    std::string workdir = workdir0;
    mkdirs(workdir);
    const std::string clean_fa = workdir + "/_" + name + ".genome.fa", clean_gff = workdir + "/_" + name + ".clean.gff";
    // 1. FASTA headers: keep the first token
    if (!nonempty(clean_fa)) {
        write_text_file(clean_fa, [&](FILE* o) {
            read_text_lines(fasta, [&](const std::string& line) {
                if (!line.empty() && line[0] == '>') {
                    std::string h = line;
                    while (!h.empty() && h.back() == '\n') h.pop_back();
                    size_t p = h.find_first_of(" \t");
                    if (p != std::string::npos) h.resize(p);
                    fputs(h.c_str(), o);
                    fputc('\n', o);
                } else fwrite(line.data(), 1, line.size(), o);
            });
        });
    }
    // 2. index + chromosome-scale sequences
    std::vector<std::pair<std::string, long long>> lens;
    const std::string fai = clean_fa + ".fai";
    if (!nonempty(fai)) {
        std::string err;
        if (!write_fai(clean_fa, fai, nullptr, &err)) die("indexing " + clean_fa + ": " + err);
    }
    std::vector<std::string> fids;
    std::unordered_map<std::string, long long> flen;
    {
        LineReader rd(fai);
        std::string_view l;
        while (rd.next(l)) {
            auto p = split_tab(l);
            if (p.size() < 2) continue;
            long long v;
            if (!py_int(p[1], v)) die(fai + ": bad length");
            std::string id(p[0]);
            if (!flen.count(id)) fids.push_back(id);
            flen[id] = v;
        }
    }
    std::unordered_set<std::string> chrom;
    for (auto& id : fids) if (flen[id] >= min_chr_len) chrom.insert(id);
    // 3. GFF seqids and their largest coordinate
    std::vector<std::string> gids;
    std::unordered_map<std::string, long long> gspan;
    {
        LineReader rd(gff);
        std::string_view l;
        while (rd.next(l)) {
            if ((!l.empty() && l[0] == '#') || py_strip(l).empty()) continue;
            auto p = split_tab(l);
            if (p.size() < 5) continue;
            long long e;
            if (!py_int(p[4], e)) continue;
            std::string sid(p[0]);
            auto it = gspan.find(sid);
            if (it == gspan.end()) { gids.push_back(sid); gspan.emplace(sid, e); }
            else if (e > it->second) it->second = e;
        }
    }
    // 4. seqid map: exact, normalised token, unique length match (within 0.5%)
    std::unordered_map<std::string, std::string> fnorm;
    for (auto& id : fids) fnorm[norm_token(id)] = id;
    std::unordered_map<std::string, std::string> smap;
    for (auto& g : gids) {
        if (flen.count(g)) { smap[g] = g; continue; }
        auto nt = fnorm.find(norm_token(g));
        if (nt != fnorm.end()) { smap[g] = nt->second; continue; }
        const long long ge = gspan[g];
        std::string bestid;
        bool amb = false;
        for (auto& id : fids) {
            long long fl = flen[id];
            if (ge <= fl && (double)(fl - ge) / (double)fl <= 0.005) {
                if (!bestid.empty()) { amb = true; break; }
                bestid = id;
            }
        }
        if (!amb && !bestid.empty()) smap[g] = bestid;
    }
    // 5. clean GFF: FASTA seqids, chromosomes only, strand "?" dropped
    write_text_file(clean_gff, [&](FILE* o) {
        LineReader rd(gff);
        std::string_view l;
        while (rd.next(l)) {
            if (!l.empty() && l[0] == '#') continue;
            auto p = split_tab(l);
            if (p.size() < 8) continue;
            auto it = smap.find(std::string(p[0]));
            if (it == smap.end() || !chrom.count(it->second)) continue;
            if (p[6] == "?") continue;
            fputs(it->second.c_str(), o);
            for (size_t i = 1; i < p.size(); ++i) { fputc('\t', o); fwrite(p[i].data(), 1, p[i].size(), o); }
            fputc('\n', o);
        }
    });
    logm(sfmt("[genome_pair:%s] chromosomes %zu/%zu (>= %s bp), seqid map %zu gff->fasta", name.c_str(), chrom.size(), fids.size(),
              commas(min_chr_len).c_str(), smap.size()));
    // 6. gffread
    const std::string log = workdir + "/_gffread_" + name + ".log";
    RunResult r = run_proc({"gffread", "-w", out_fa, "-g", clean_fa, clean_gff}, log, "&1");
    if (r.status != 0) die("gffread failed (exit " + std::to_string(r.status) + "), see " + log);
    // 7. ">PREFIX__id"
    if (!id_prefix.empty()) {
        const std::string tmp = out_fa + ".tmp";
        FILE* o = fopen(tmp.c_str(), "w");
        if (!o) die("cannot write " + tmp);
        read_text_lines(out_fa, [&](const std::string& line) {
            if (!line.empty() && line[0] == '>') {
                auto w = py_split_ws(std::string_view(line).substr(1));
                if (w.empty()) die(out_fa + ": FASTA header without an id");
                fprintf(o, ">%s__%s\n", id_prefix.c_str(), std::string(w[0]).c_str());
            } else fwrite(line.data(), 1, line.size(), o);
        });
        if (fclose(o) != 0 || rename(tmp.c_str(), out_fa.c_str()) != 0) die("cannot write " + out_fa);
    }
}

int extract_tx_main(int argc, char** argv) {
    std::string fasta, gff, workdir, name = "genome", out, prefix;
    long long minlen = 1000000;
    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&]() -> std::string {
            if (i + 1 >= argc) die("missing value for " + a);
            return argv[++i];
        };
        if (a == "--fasta") fasta = need();
        else if (a == "--gff") gff = need();
        else if (a == "--workdir") workdir = need();
        else if (a == "--name") name = need();
        else if (a == "--out") out = need();
        else if (a == "--id-prefix") prefix = need();
        else if (a == "--min-chr-len") { if (!py_int(need(), minlen)) die("--min-chr-len: not an integer"); }
        else if (a == "-h" || a == "--help") {
            std::fprintf(stderr,
                         "usage: lq extract-transcripts --fasta GENOME.fa --gff GENOME.gff --workdir DIR --out TX.fa\n"
                         "                                        [--name genome] [--id-prefix P] [--min-chr-len 1000000]\n"
                         "  FASTA headers cut at the first space, sequences >= --min-chr-len kept, GFF seqids matched to the\n"
                         "  FASTA (exact, chr01 ~ Chr1 ~ CM0001.1 by name, or by length within 0.5%%), strand '?' dropped,\n"
                         "  then gffread -w. --id-prefix P writes >P__<id>. Same as utils/genome_pair.py.\n");
            return 0;
        } else die("unknown argument " + a);
    }
    if (fasta.empty() || gff.empty() || workdir.empty() || out.empty()) die("--fasta, --gff, --workdir and --out are required");
    genome_pair_extract(fasta, gff, workdir, name, out, prefix, minlen);
    long long n = 0;
    LineReader rd(out);
    std::string_view l;
    while (rd.next(l)) if (!l.empty() && l[0] == '>') ++n;
    std::printf("[done] %lld transcripts -> %s\n", n, out.c_str());
    return 0;
}

// ===================================================================================== normalize-clusters
namespace {
std::string clean_id(std::string_view g) { return std::string(before_dot(py_strip(g))); }
struct UF {   // union-find; components in the order of their first-inserted member, members in insertion order
    std::unordered_map<std::string, uint32_t> ix;
    std::vector<std::string> nodes;
    std::vector<uint32_t> par;
    uint32_t id(const std::string& s) {
        auto it = ix.find(s);
        if (it != ix.end()) return it->second;
        uint32_t i = (uint32_t)nodes.size();
        ix.emplace(s, i);
        nodes.push_back(s);
        par.push_back(i);
        return i;
    }
    uint32_t find(uint32_t x) {
        while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; }
        return x;
    }
    void unite(const std::string& a, const std::string& b) {
        uint32_t ra = find(id(a)), rb = find(id(b));
        if (ra != rb) par[rb] = ra;
    }
    std::vector<std::vector<std::string>> components() {
        std::unordered_map<uint32_t, uint32_t> cix;
        std::vector<std::vector<std::string>> out;
        for (uint32_t i = 0; i < nodes.size(); ++i) {
            uint32_t r = find(i);
            auto it = cix.find(r);
            if (it == cix.end()) { cix.emplace(r, (uint32_t)out.size()); out.push_back({nodes[i]}); }
            else out[it->second].push_back(nodes[i]);
        }
        return out;
    }
};
}  // namespace

int normalize_clusters_main(int argc, char** argv) {
    std::string input, format, clusters_out, ancestry_out, orphan_out, all_genes, id_col;
    std::vector<std::string> tags, tsv_cols;
    bool has_id_col = false;
    int min_size = 2;
    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&]() -> std::string {
            if (i + 1 >= argc) die("missing value for " + a);
            return argv[++i];
        };
        if (a == "--input") input = need();
        else if (a == "--format") format = need();
        else if (a == "--clusters-out") clusters_out = need();
        else if (a == "--ancestry-out") ancestry_out = need();
        else if (a == "--orphan-out") orphan_out = need();
        else if (a == "--all-genes") all_genes = need();
        else if (a == "--subgenome-tag") tags.push_back(need());
        else if (a == "--tsv-cols") { while (i + 1 < argc && argv[i + 1][0] != '-') tsv_cols.push_back(argv[++i]); }
        else if (a == "--tsv-id-col") { id_col = need(); has_id_col = true; }
        else if (a == "--min-size") min_size = std::atoi(need().c_str());
        else if (a == "-h" || a == "--help") {
            std::fprintf(stderr,
                         "usage: lq normalize-clusters --input RAW --format standard|jcvi|tsv|pairs --clusters-out OUT\n"
                         "        [--ancestry-out FILE --subgenome-tag NAME:REGEX ...] [--orphan-out FILE --all-genes FILE]\n"
                         "        [--tsv-cols COL ...] [--tsv-id-col COL] [--min-size 2]\n"
                         "  standard: ClusterID<TAB>g1,g2,...   jcvi / pairs: gene pairs -> connected components\n"
                         "  tsv: one cluster per row across columns (names or 0-based indices). Gene ids are cut at the first '.'.\n"
                         "  Output: Global_Cluster_NNNNNN<TAB>sorted,unique,genes (>= --min-size). Same as normalize_clusters.py\n"
                         "  (regexes: ECMAScript syntax, which covers the usual ^PREFIX patterns).\n");
            return 0;
        } else die("unknown argument " + a);
    }
    if (input.empty() || format.empty() || clusters_out.empty()) die("--input, --format and --clusters-out are required");
    std::vector<std::vector<std::string>> clusters;
    if (format == "standard") {
        LineReader rd(input);
        std::string_view l;
        while (rd.next(l)) {
            if (py_strip(l).empty() || (!l.empty() && l[0] == '#')) continue;
            auto p = py_split(l, '\t');
            if (p.size() < 2) continue;
            std::vector<std::string> g;
            for (auto x : py_split(p[1], ',')) if (!py_strip(x).empty()) g.push_back(clean_id(x));
            if (!g.empty()) clusters.push_back(std::move(g));
        }
    } else if (format == "jcvi" || format == "pairs") {
        UF uf;
        long long n = 0;
        LineReader rd(input);
        std::string_view l;
        while (rd.next(l)) {
            if (py_strip(l).empty() || (!l.empty() && l[0] == '#')) continue;
            auto p = py_split(l, '\t');
            if (p.size() < 2) p = py_split_ws(l);
            if (p.size() < 2) continue;
            std::string a = clean_id(p[0]), b = clean_id(p[1]);
            if (!a.empty() && !b.empty()) { uf.unite(a, b); ++n; }
        }
        logm("[normalize] parsed " + std::to_string(n) + " pairs -> components");
        clusters = uf.components();
    } else if (format == "tsv") {
        LineReader rd(input);
        std::string_view l;
        if (!rd.next(l)) l = std::string_view();
        std::vector<std::string> header;
        for (auto x : py_split(l, '\t')) header.emplace_back(x);
        auto col_of = [&](const std::string& c) -> long long {
            for (size_t i = 0; i < header.size(); ++i) if (header[i] == c) return (long long)i;
            long long v;
            if (!py_int(c, v)) die("--tsv-cols / --tsv-id-col: '" + c + "' is neither a column name nor an index");
            return v;
        };
        std::vector<long long> cols;
        if (tsv_cols.empty()) for (size_t i = 0; i < header.size(); ++i) cols.push_back((long long)i);
        else for (auto& c : tsv_cols) cols.push_back(col_of(c));
        while (rd.next(l)) {
            if (py_strip(l).empty()) continue;
            auto p = py_split(l, '\t');
            std::vector<std::string> g;
            for (long long ci : cols) {
                long long cj = ci < 0 ? ci + (long long)p.size() : ci;   // Python negative index
                if (ci < (long long)p.size() && cj >= 0 && !py_strip(p[cj]).empty() && p[cj] != ".") g.push_back(clean_id(p[cj]));
            }
            if (has_id_col) {   // as normalize_clusters.py: genes[j] is compared with the j-th requested column
                long long idc = col_of(id_col);
                std::vector<std::string> g2;
                for (size_t j = 0; j < g.size(); ++j) if (cols[j] != idc) g2.push_back(g[j]);
                g.swap(g2);
            }
            if (!g.empty()) clusters.push_back(std::move(g));
        }
    } else die("--format must be standard, jcvi, tsv or pairs");

    // clusters
    long long n = 0;
    write_text_file(clusters_out, [&](FILE* f) {
        for (auto& c : clusters) {
            std::set<std::string> u(c.begin(), c.end());
            if ((int)u.size() < min_size) continue;
            std::string j;
            for (auto& g : u) j += (j.empty() ? "" : ",") + g;
            fprintf(f, "Global_Cluster_%06lld\t%s\n", n, j.c_str());
            ++n;
        }
    });
    logm("[normalize] wrote " + std::to_string(n) + " clusters (>= " + std::to_string(min_size) + ") -> " + clusters_out);
    // ancestry
    if (!ancestry_out.empty()) {
        std::vector<std::pair<std::string, std::regex>> rx;
        for (auto& t : tags) {
            size_t c = t.find(':');
            if (c == std::string::npos) die("--subgenome-tag must be NAME:REGEX, got " + t);
            try { rx.push_back({t.substr(0, c), std::regex(t.substr(c + 1), std::regex::ECMAScript)}); }
            catch (const std::regex_error& e) { die("--subgenome-tag " + t + ": bad regex (" + e.what() + ")"); }
        }
        if (rx.empty()) logm("[normalize] no --subgenome-tag; all lineages UNKNOWN");
        std::set<std::string> genes;
        for (auto& c : clusters) {
            std::set<std::string> u(c.begin(), c.end());
            if (u.size() >= 2) for (auto& g : c) genes.insert(clean_id(g));
        }
        long long unk = 0;
        write_text_file(ancestry_out, [&](FILE* f) {
            fputs("GeneID\tAncestry\n", f);
            for (auto& g : genes) {
                std::string lin = "UNKNOWN";
                for (auto& r : rx) if (std::regex_search(g, r.second)) { lin = r.first; break; }
                if (lin == "UNKNOWN") ++unk;
                fprintf(f, "%s\t%s\n", g.c_str(), lin.c_str());
            }
        });
        if (unk) logm("[normalize] WARNING: " + std::to_string(unk) + " genes unmatched by any --subgenome-tag (lineage=UNKNOWN)");
        logm("[normalize] wrote ancestry for " + std::to_string(genes.size()) + " genes -> " + ancestry_out);
    }
    // orphans
    if (!orphan_out.empty()) {
        std::set<std::string> in_cluster, all;
        for (auto& c : clusters) {
            std::set<std::string> u(c.begin(), c.end());
            if (u.size() >= 2) for (auto& g : c) in_cluster.insert(clean_id(g));
        }
        if (!all_genes.empty()) {
            LineReader rd(all_genes);
            std::string_view l;
            while (rd.next(l)) if (!py_strip(l).empty()) all.insert(clean_id(l));
        }
        long long i = 0;
        write_text_file(orphan_out, [&](FILE* f) {
            for (auto& g : all) if (!in_cluster.count(g)) fprintf(f, "OrphanGrp_%06lld\t%s\n", ++i, g.c_str());
        });
        logm("[normalize] wrote " + std::to_string(i) + " orphan genes -> " + orphan_out);
    }
    return 0;
}

}  // namespace lq
