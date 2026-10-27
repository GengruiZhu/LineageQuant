// lqcore -- minimal SAM record parsing (no allocation)
#pragma once
#include "common.hpp"

namespace lq {

struct Rec {
    std::string_view qname, rname, cigar, seq;
    int flag = 0, nm = 0;
    long long pos = 0;      // 1-based POS as in SAM
    long long pnext = 0;    // 1-based PNEXT
    long long tl = 0;       // |TLEN|
    uint32_t rid = 0xFFFFFFFFu;
};

// QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAGS]
inline bool parse_sam(std::string_view l, Rec& r) {
    size_t st[11];
    st[0] = 0;
    size_t p = 0;
    int nf = 1;
    while (nf < 11) {
        size_t q = l.find('\t', p);
        if (q == std::string_view::npos) return false;
        st[nf++] = q + 1;
        p = q + 1;
    }
    auto fld = [&](int i) -> std::string_view {
        size_t b = st[i];
        size_t e = l.find('\t', b);
        if (e == std::string_view::npos) e = l.size();
        return l.substr(b, e - b);
    };
    r.qname = l.substr(0, st[1] - 1);
    r.flag = (int)parse_ll(fld(1));
    r.rname = fld(2);
    r.pos = parse_ll(fld(3));
    r.cigar = fld(5);
    r.pnext = parse_ll(fld(7));
    r.tl = llabs(parse_ll(fld(8)));
    r.seq = fld(9);
    r.nm = 0;
    size_t e = l.find('\t', st[10]);   // end of QUAL -> tags
    while (e != std::string_view::npos) {
        size_t nb = e + 1;
        size_t ne = l.find('\t', nb);
        std::string_view tg = l.substr(nb, (ne == std::string_view::npos ? l.size() : ne) - nb);
        if (tg.size() > 5 && tg[0] == 'N' && tg[1] == 'M' && tg[2] == ':' && tg[3] == 'i' && tg[4] == ':') {
            r.nm = (int)parse_ll(tg.substr(5));
            break;
        }
        e = ne;
    }
    return true;
}

// iterate CIGAR operations: f(op, len)
template <class F>
inline void for_cigar(std::string_view cg, F&& f) {
    long long n = 0;
    for (char c : cg) {
        if (c >= '0' && c <= '9') n = n * 10 + (c - '0');
        else { f(c, n); n = 0; }
    }
}

}  // namespace lq
