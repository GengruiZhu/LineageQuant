// LineageQuant 1.1.0 -- binary outputs
//
//   Theta    (theta --binary)      DIR/Global_Theta_Matrix.lqb      magic "LQTHETA1"
//   result   (quant -o S.lqres)    EM + effective lengths + per-reference counts of one sample, magic "LQRES001"
//
// Both hold exactly what the text files hold (the same numbers before printing); `lq export` writes the
// text files byte for byte as the text writers do, and every reader that takes the text takes the binary too.
// Little-endian, fixed-width integers; strings = u32 length + bytes. Layout: docs/FORMATS.md.
#pragma once
#include "common.hpp"

#include <string>
#include <vector>

namespace lq {

// ------------------------------------------------------------------ byte buffers
struct BinOut {
    std::vector<char> b;
    void raw(const void* p, size_t n) { b.insert(b.end(), (const char*)p, (const char*)p + n); }
    void u8(uint8_t v) { raw(&v, 1); }
    void u32(uint32_t v) { raw(&v, 4); }
    void u64(uint64_t v) { raw(&v, 8); }
    void i64(int64_t v) { raw(&v, 8); }
    void f64(double v) { raw(&v, 8); }
    void str(const std::string& s) { u32((uint32_t)s.size()); raw(s.data(), s.size()); }
    // tmp file + rename: a reader never sees a half-written file
    void save(const std::string& path) const {
        std::string tmp = path + ".tmp" + std::to_string((int)getpid());
        FILE* f = fopen(tmp.c_str(), "wb");
        if (!f) die("cannot write " + tmp);
        if ((!b.empty() && fwrite(b.data(), 1, b.size(), f) != b.size()) || fclose(f) != 0) {
            unlink(tmp.c_str());
            die("error writing " + path);
        }
        if (rename(tmp.c_str(), path.c_str()) != 0) die("cannot rename " + tmp + " -> " + path);
    }
};
struct BinIn {
    std::vector<char> b;
    size_t p = 0;
    std::string path;
    explicit BinIn(const std::string& f) : path(f) {
        FILE* fp = fopen(f.c_str(), "rb");
        if (!fp) die("cannot open " + f);
        char buf[1 << 16];
        size_t n;
        while ((n = fread(buf, 1, sizeof buf, fp)) > 0) b.insert(b.end(), buf, buf + n);
        fclose(fp);
    }
    void need(size_t n) const { if (p + n > b.size()) die(path + ": truncated binary file"); }
    void raw(void* d, size_t n) { need(n); if (n) memcpy(d, b.data() + p, n); p += n; }
    uint8_t u8() { uint8_t v; raw(&v, 1); return v; }
    uint32_t u32() { uint32_t v; raw(&v, 4); return v; }
    uint64_t u64() { uint64_t v; raw(&v, 8); return v; }
    int64_t i64() { int64_t v; raw(&v, 8); return v; }
    double f64() { double v; raw(&v, 8); return v; }
    std::string str() { uint32_t n = u32(); need(n); std::string s(b.data() + p, n); p += n; return s; }
};
inline bool has_magic(const std::string& path, const char* magic8) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    char m[8];
    bool ok = fread(m, 1, 8, f) == 8 && memcmp(m, magic8, 8) == 0;
    fclose(f);
    return ok;
}

// ------------------------------------------------------------------ Theta
static const char THETA_MAGIC[9] = "LQTHETA1";
struct ThetaRow {
    uint32_t cl = 0;                 // index into ThetaBin::clusters
    std::string allele;
    uint64_t probes = 0, hits = 0;
    uint8_t flags = 0;               // THF_*
    double lr = 0.0;                 // long-read count (if THF_LR)
    double theta = 0.0;              // unrounded; the text holds "%.6f"
    double shrink = 0.0;             // shrinkage used (Kmer_Inferred)
};
enum : uint8_t { THF_LR = 1, THF_ANCHORED = 2 };
struct ThetaBin {
    std::vector<std::string> clusters;
    std::vector<ThetaRow> rows;
};
inline void write_theta_bin(const std::string& path, const ThetaBin& T) {
    BinOut o;
    o.raw(THETA_MAGIC, 8);
    o.u32((uint32_t)T.clusters.size());
    for (auto& c : T.clusters) o.str(c);
    o.u64(T.rows.size());
    for (auto& r : T.rows) {
        o.u32(r.cl); o.str(r.allele); o.u64(r.probes); o.u64(r.hits); o.u8(r.flags);
        o.f64(r.lr); o.f64(r.theta); o.f64(r.shrink);
    }
    o.save(path);
}
inline void read_theta_bin(const std::string& path, ThetaBin& T) {
    BinIn in(path);
    char m[8];
    in.raw(m, 8);
    if (memcmp(m, THETA_MAGIC, 8) != 0) die(path + ": not a LineageQuant Theta file");
    T.clusters.resize(in.u32());
    for (auto& c : T.clusters) c = in.str();
    T.rows.resize(in.u64());
    for (auto& r : T.rows) {
        r.cl = in.u32(); r.allele = in.str(); r.probes = in.u64(); r.hits = in.u64(); r.flags = in.u8();
        r.lr = in.f64(); r.theta = in.f64(); r.shrink = in.f64();
        if (r.cl >= T.clusters.size()) die(path + ": corrupt (cluster index)");
    }
}
// the Long-read column: "0" without a long-read count, else Python repr of the count
inline std::string theta_lr_str(const ThetaRow& r) { return (r.flags & THF_LR) ? py_repr(r.lr) : std::string("0"); }
inline void theta_tsv_header(FILE* f) {
    std::fputs("Cluster_ID\tAllele_ID\tUnique_Probes\tRaw_Hits\tLR_Count\tTheta\tPrior_Status\tShrinkage_Used\n", f);
}
inline void theta_tsv_row(FILE* f, const ThetaBin& T, const ThetaRow& r) {
    const std::string ls = theta_lr_str(r);
    if (r.flags & THF_ANCHORED)
        std::fprintf(f, "%s\t%s\t%llu\t%llu\t%s\t%.6f\tLongRead_Anchored\tNA\n", T.clusters[r.cl].c_str(), r.allele.c_str(),
                     (unsigned long long)r.probes, (unsigned long long)r.hits, ls.c_str(), r.theta);
    else
        std::fprintf(f, "%s\t%s\t%llu\t%llu\t%s\t%.6f\tKmer_Inferred\t%.2f\n", T.clusters[r.cl].c_str(), r.allele.c_str(),
                     (unsigned long long)r.probes, (unsigned long long)r.hits, ls.c_str(), r.theta, r.shrink);
}

// ------------------------------------------------------------------ result of one sample
static const char RES_MAGIC[9] = "LQRES001";
struct ResultBin {
    int32_t iter = 0;
    // targets (EM rows, in EM order)
    std::vector<std::string> cluster, allele, ancestry;
    std::vector<double> th0, counts, fin;    // unrounded; the text holds %.6f / %.2f / %.6f
    std::vector<double> tlen, eff, gc;       // effective lengths the EM used (text: %.0f / %.4f / %.5f)
    // references (samtools idxstats semantics)
    std::vector<std::string> refs;
    std::vector<int64_t> rlen;
    std::vector<uint64_t> mapped, unmapped;
    uint64_t unmapped_star = 0;
    size_t T() const { return allele.size(); }
};
inline void write_result_bin(const std::string& path, const ResultBin& R) {
    BinOut o;
    o.b.reserve(64 + R.T() * 96 + R.refs.size() * 48);
    o.raw(RES_MAGIC, 8);
    o.i64(R.iter);
    o.u64(R.T());
    for (size_t i = 0; i < R.T(); ++i) {
        o.str(R.cluster[i]); o.str(R.allele[i]); o.str(R.ancestry[i]);
        o.f64(R.th0[i]); o.f64(R.counts[i]); o.f64(R.fin[i]); o.f64(R.tlen[i]); o.f64(R.eff[i]); o.f64(R.gc[i]);
    }
    o.u64(R.refs.size());
    for (size_t i = 0; i < R.refs.size(); ++i) { o.str(R.refs[i]); o.i64(R.rlen[i]); o.u64(R.mapped[i]); o.u64(R.unmapped[i]); }
    o.u64(R.unmapped_star);
    o.save(path);
}
inline void read_result_bin(const std::string& path, ResultBin& R) {
    BinIn in(path);
    char m[8];
    in.raw(m, 8);
    if (memcmp(m, RES_MAGIC, 8) != 0) die(path + ": not a LineageQuant result file");
    R.iter = (int32_t)in.i64();
    size_t T = in.u64();
    for (auto* v : {&R.cluster, &R.allele, &R.ancestry}) v->resize(T);
    for (auto* v : {&R.th0, &R.counts, &R.fin, &R.tlen, &R.eff, &R.gc}) v->resize(T);
    for (size_t i = 0; i < T; ++i) {
        R.cluster[i] = in.str(); R.allele[i] = in.str(); R.ancestry[i] = in.str();
        R.th0[i] = in.f64(); R.counts[i] = in.f64(); R.fin[i] = in.f64(); R.tlen[i] = in.f64(); R.eff[i] = in.f64(); R.gc[i] = in.f64();
    }
    size_t NR = in.u64();
    R.refs.resize(NR); R.rlen.resize(NR); R.mapped.resize(NR); R.unmapped.resize(NR);
    for (size_t i = 0; i < NR; ++i) { R.refs[i] = in.str(); R.rlen[i] = in.i64(); R.mapped[i] = in.u64(); R.unmapped[i] = in.u64(); }
    R.unmapped_star = in.u64();
}
// text writers (the only place these formats are printed)
inline void write_em_tsv(FILE* f, const ResultBin& R) {
    fprintf(f, "ClusterID\tAlleleID\tAncestry\tInitial_Theta\tEM_Allocated_Reads\tFinal_Theta\tEM_Iter\tEM_LogLik\n");
    for (size_t i = 0; i < R.T(); ++i)
        fprintf(f, "%s\t%s\t%s\t%.6f\t%.2f\t%.6f\t%d\t0\n", R.cluster[i].c_str(), R.allele[i].c_str(), R.ancestry[i].c_str(),
                R.th0[i], R.counts[i], R.fin[i], (int)R.iter);
}
inline void write_eff_tsv(FILE* f, const ResultBin& R) {
    fprintf(f, "Target\tLength\tEffectiveLength\tGC_Factor\n");
    for (size_t i = 0; i < R.T(); ++i) fprintf(f, "%s\t%.0f\t%.4f\t%.5f\n", R.allele[i].c_str(), R.tlen[i], R.eff[i], R.gc[i]);
}
inline void write_idx_tsv(FILE* f, const ResultBin& R) {
    for (size_t i = 0; i < R.refs.size(); ++i)
        fprintf(f, "%s\t%lld\t%llu\t%llu\n", R.refs[i].c_str(), (long long)R.rlen[i], (unsigned long long)R.mapped[i],
                (unsigned long long)R.unmapped[i]);
    fprintf(f, "*\t0\t0\t%llu\n", (unsigned long long)R.unmapped_star);
}
// writes a text file through tmp + rename
template <class F>
inline void write_text_file(const std::string& path, F&& body) {
    std::string tmp = path + ".tmp" + std::to_string((int)getpid());
    FILE* f = fopen(tmp.c_str(), "w");
    if (!f) die("cannot write " + tmp);
    std::vector<char> obuf(1 << 22);
    setvbuf(f, obuf.data(), _IOFBF, obuf.size());
    body(f);
    if (fclose(f) != 0) { unlink(tmp.c_str()); die("error writing " + path); }
    if (rename(tmp.c_str(), path.c_str()) != 0) die("cannot rename " + tmp + " -> " + path);
}
inline bool ends_with(const std::string& s, const std::string& x) {
    return s.size() >= x.size() && s.compare(s.size() - x.size(), x.size(), x) == 0;
}

}  // namespace lq
