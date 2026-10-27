// lqcore -- multithreaded, order-preserving scan of a name-grouped SAM stream
#include "scan.hpp"

#include <cmath>
#include <condition_variable>
#include <deque>
#include <map>
#include <mutex>
#include <thread>

namespace lq {

// ================================================================== P table (em -p)
void PTable::load(const std::string& path, const Opts& o) {
    LineReader rd(path);
    std::string_view l;
    if (!rd.next(l)) l = std::string_view();
    auto cols = split_tab(l);
    std::vector<size_t> pci;
    int si = -1;
    for (size_t i = 0; i < cols.size(); ++i) {
        std::string c(cols[i]);
        size_t b = 0, e = c.size();
        while (b < e && isspace((unsigned char)c[b])) ++b;
        while (e > b && isspace((unsigned char)c[e - 1])) --e;
        std::string cs = c.substr(b, e - b);
        if (cs.size() >= 4 && cs[0] == 'P' && cs[1] == '(' && cs.back() == ')') {
            lineages.push_back(upper(cs.substr(2, cs.size() - 3)));
            pci.push_back(i);
        }
        std::string lc = c;
        for (auto& ch : lc) ch = (char)tolower((unsigned char)ch);
        if (lc.find("status") != std::string::npos) si = (int)i;
    }
    L = pci.size();
    const double WC[3] = {o.pw_spec, o.pw_conf, o.pw_none};
    std::vector<double> pv(L);
    auto cls = [](std::string_view s) {
        auto ist = [&](std::string_view p) {
            if (s.size() < p.size()) return false;
            for (size_t i = 0; i < p.size(); ++i) if (tolower((unsigned char)s[i]) != p[i]) return false;
            return true;
        };
        if (ist("specific")) return 0;
        std::string lc(s);
        for (auto& c : lc) c = (char)tolower((unsigned char)c);
        for (const char* k : {"none", "no_snp", "chimeric", "ambiguous", "weak"})
            if (lc.find(k) != std::string::npos) return 2;
        return 1;
    };
    while (rd.next(l)) {
        auto q = split_tab(l);
        std::string rn(q[0]);
        if (idx.count(rn)) continue;
        int wc = (si >= 0 && (int)q.size() > si) ? cls(q[si]) : 1;
        double w = WC[wc];
        if (w == 0.0) continue;
        bool ok = true;
        for (size_t j = 0; j < L && ok; ++j)
            if (pci[j] >= q.size() || !py_float(q[pci[j]], pv[j])) ok = false;
        if (!ok) continue;
        for (size_t j = 0; j < L; ++j) {
            double x = pv[j];
            if (o.p_alpha > 0) x = (x + o.p_alpha) / (1.0 + (double)L * o.p_alpha);
            x = std::max(x, o.p_floor);
            if (w != 1.0) x = std::pow(x, w);
            pv[j] = x;
        }
        idx.emplace(rn, (uint32_t)(val.size() / (L ? L : 1)));
        val.insert(val.end(), pv.begin(), pv.end());
    }
    std::string ls = "[";
    for (size_t i = 0; i < lineages.size(); ++i) ls += (i ? ", '" : "'") + lineages[i] + "'";
    ls += "]";
    logm("[INFO] lineage P: " + ls + "; reads with a non-zero weight: " + commas((long long)idx.size()));
}

// ================================================================== header
bool read_sam_header(LineReader& rd, std::vector<std::string>& refs, std::vector<long long>& lens, std::string& first) {
    std::string_view line;
    while (rd.next(line)) {
        if (line.rfind("@", 0) == 0) {
            if (line.rfind("@SQ", 0) == 0) {
                std::string sn;
                long long ln = -1;
                auto f = split_tab(line);
                for (size_t i = 1; i < f.size(); ++i) {
                    size_t c = f[i].find(':');
                    if (c == std::string_view::npos) continue;
                    auto key = f[i].substr(0, c), v = f[i].substr(c + 1);
                    if (key == "SN") sn = std::string(v);
                    else if (key == "LN") ln = parse_ll(v);
                }
                refs.push_back(sn);
                lens.push_back(ln);
            }
            continue;
        }
        if (line.empty()) continue;
        first.assign(line);
        return true;
    }
    return false;
}

// ================================================================== scan
namespace {
using Batch = ScanBatch;
struct BatchOut {
    std::vector<uint8_t> kind;  // 1 single, 2 eq, 3 eqx
    std::vector<uint32_t> tgt;
    std::string keys;
    std::vector<uint64_t> koff;
    long long n_rec = 0, n_frag = 0, n_unal = 0, nm_sum = 0, len_sum = 0;
    std::vector<int32_t> fl;
    std::vector<uint32_t> gcb;  // 1.1.0: fragments per GC percent
    std::string dk, ds, pr, df;
};
struct Acc {
    std::vector<uint64_t> im, iu;
    uint64_t star = 0;
    LineageStats ls;
};
struct Ent { uint32_t t; int e0, e1; };
}  // namespace

void scan_sam(LineReader& rd, const std::string& first, const ScanCfg& cfg, Evidence& ev, const ScanPrepare& prepare) {
    const size_t BATCH_BYTES = 8u << 20;
    ScanProducer produce = [&](const ScanEmit& emit) {
        std::unique_ptr<Batch> b(new Batch);
        std::string curq;
        bool have = false, keep = true;
        long long fidx = -1;
        auto flush = [&]() {
            if (b->fs.empty()) return;
            b->ls.push_back((uint32_t)b->data.size());
            b->fs.push_back((uint32_t)(b->ls.size() - 1));
            emit(std::move(b));
            b.reset(new Batch);
        };
        auto add = [&](std::string_view l) {
            size_t t = l.find('\t');
            std::string_view q = l.substr(0, t);
            if (!have || q != curq) {
                if (b->data.size() >= BATCH_BYTES) flush();
                curq.assign(q);
                have = true;
                ++fidx;
                keep = cfg.chunk_n <= 1 || ((fidx / cfg.block) % cfg.chunk_n == cfg.chunk_i - 1);
                if (keep) b->fs.push_back((uint32_t)b->ls.size());
            }
            if (keep) {
                if (b->data.size() + l.size() + 1 >= 0xFFFFFFF0u) die("batch overflow");
                b->ls.push_back((uint32_t)b->data.size());
                b->data.append(l);
                b->data.push_back('\n');
            }
        };
        if (!first.empty()) add(first);
        std::string_view l;
        while (rd.next(l))
            if (!l.empty()) add(l);
        flush();
    };
    scan_batches(produce, prepare, cfg, ev);
}

void scan_batches(const ScanProducer& produce, const ScanPrepare& prepare, const ScanCfg& cfg, Evidence& ev) {
    const Model& M = *cfg.M;
    const Opts& o = *cfg.O;
    const bool NMW = o.nm_weight != "off";
    const int W = std::max(1, cfg.workers);
    const size_t MAX_INFLIGHT = 4 * (size_t)W + 2;
    const size_t L = cfg.LE ? cfg.LE->em_lins.size() : (cfg.PT ? cfg.PT->L : 0);

    std::mutex mu;
    std::condition_variable cv_work, cv_res, cv_space;
    std::deque<std::unique_ptr<Batch>> work;
    std::map<size_t, std::unique_ptr<BatchOut>> results;
    bool reader_done = false;
    size_t produced = 0, merged = 0;

    // ---------------------------------------------------------------- reader
    std::thread reader([&]() {
        ScanEmit emit = [&](std::unique_ptr<Batch> b) {
            {
                std::unique_lock<std::mutex> lk(mu);
                cv_space.wait(lk, [&]() { return produced - merged < MAX_INFLIGHT; });
                b->id = produced++;
                work.push_back(std::move(b));
            }
            cv_work.notify_one();
        };
        produce(emit);
        {
            std::lock_guard<std::mutex> lk(mu);
            reader_done = true;
        }
        cv_work.notify_all();
        cv_res.notify_all();
    });

    // ---------------------------------------------------------------- workers
    std::vector<Acc> accs(W);
    for (auto& a : accs) { a.im.assign(M.refs.size(), 0); a.iu.assign(M.refs.size(), 0); }
    auto worker = [&](int wid) {
        Acc& acc = accs[wid];
        LinScratch lsc;
        std::vector<Rec> recs;
        std::vector<Ent> ents;
        std::vector<int> tot;
        std::vector<double> mbuf;
        std::string key, rbuf, last_rname;
        uint32_t last_rid = NOCL;
        auto rid_of = [&](std::string_view rn) -> uint32_t {
            if (last_rid != NOCL && rn == last_rname) return last_rid;
            rbuf.assign(rn);
            auto it = M.r2i.find(rbuf);
            if (it == M.r2i.end()) die("reference not in header: " + rbuf);
            last_rname.assign(rn);
            last_rid = it->second;
            return last_rid;
        };
        for (;;) {
            std::unique_ptr<Batch> b;
            {
                std::unique_lock<std::mutex> lk(mu);
                cv_work.wait(lk, [&]() { return !work.empty() || reader_done; });
                if (work.empty()) return;
                b = std::move(work.front());
                work.pop_front();
            }
            if (prepare) prepare(*b, wid);
            std::unique_ptr<BatchOut> out(new BatchOut);
            const size_t nf = b->fs.empty() ? 0 : b->fs.size() - 1;
            for (size_t f = 0; f < nf; ++f) {
                recs.clear();
                for (uint32_t li = b->fs[f]; li < b->fs[f + 1]; ++li) {
                    std::string_view l(b->data.data() + b->ls[li], b->ls[li + 1] - b->ls[li] - 1);
                    Rec r;
                    if (!parse_sam(l, r)) die("malformed SAM line (fewer than 11 fields): " + std::string(l.substr(0, 200)));
                    r.rid = (r.flag & 4) ? NOCL : rid_of(r.rname);
                    recs.push_back(r);
                }
                out->n_rec += (long long)recs.size();
                out->n_frag += 1;
                long long qlen = 0;
                ents.clear();
                for (const Rec& r : recs) {
                    const bool prim = !(r.flag & 0x900);
                    if (r.seq != "*" && prim) qlen += (long long)r.seq.size();
                    if (r.rid != NOCL) {
                        acc.im[r.rid]++;
                        uint32_t t = M.ref2t[r.rid];
                        Ent* e = nullptr;
                        if (!ents.empty() && ents.back().t == t) e = &ents.back();
                        else {
                            for (auto& x : ents) if (x.t == t) { e = &x; break; }
                            if (!e) { ents.push_back({t, NM_INF, NM_INF}); e = &ents.back(); }
                        }
                        if (NMW) {
                            if (r.flag & 0x80) e->e1 = std::min(e->e1, r.nm);
                            else e->e0 = std::min(e->e0, r.nm);
                        }
                    } else if (r.rname == "*" || r.rname.empty()) acc.star++;
                    else acc.iu[rid_of(r.rname)]++;
                    if ((r.flag & 2) && prim && r.tl > 0 && r.tl < 2000) out->fl.push_back((int32_t)r.tl);
                }
                if (ents.empty()) { out->n_unal++; continue; }
                if (cfg.gc) {   // GC of the fragment = reference interval of the primary proper pair (mate 1 record)
                    if (out->gcb.empty()) out->gcb.assign(101, 0);
                    for (const Rec& r : recs) {
                        if ((r.flag & 0x40) && (r.flag & 2) && !(r.flag & 0x904) && r.rid != NOCL && r.tl > 0 && r.tl < 2000) {
                            const long long s0 = std::min(r.pos, r.pnext > 0 ? r.pnext : r.pos) - 1;
                            if (s0 >= 0 && cfg.gc->ok(r.rid) && s0 + r.tl <= cfg.gc->len[r.rid]) {
                                const uint64_t g = cfg.gc->count(r.rid, (uint64_t)s0, (uint64_t)r.tl);
                                out->gcb[(200 * g + (uint64_t)r.tl) / (2 * (uint64_t)r.tl)]++;
                            }
                            break;
                        }
                    }
                }
                int best = 0;
                if (NMW) {
                    best = INT32_MAX;
                    for (auto& e : ents) best = std::min(best, (e.e0 < NM_INF ? e.e0 : 0) + (e.e1 < NM_INF ? e.e1 : 0));
                    if (qlen) { out->nm_sum += best; out->len_sum += qlen; }
                }
                // lineage P for this fragment
                double pvb[8];
                const double* pv = nullptr;
                if (cfg.LE) {
                    if (cfg.LE->eval(recs.data(), recs.size(), recs[0].qname, lsc, acc.ls, pvb,
                                     cfg.dbg_kmer ? &out->dk : nullptr, cfg.dbg_snp ? &out->ds : nullptr,
                                     cfg.pvals_out ? &out->pr : nullptr))
                        pv = pvb;
                }
                if (cfg.dbg_frag) {          // Read_Name  N  target:extra_mismatches,...  P(lineages) | NA
                    std::string& d = out->df;
                    d.append(recs[0].qname);
                    d.push_back('\t');
                    d.append(std::to_string(ents.size()));
                    d.push_back('\t');
                    for (size_t i = 0; i < ents.size(); ++i) {
                        if (i) d.push_back(',');
                        d.append(M.names[ents[i].t]);
                        d.push_back(':');
                        int dn = NMW ? (ents[i].e0 < NM_INF ? ents[i].e0 : 0) + (ents[i].e1 < NM_INF ? ents[i].e1 : 0) - best : 0;
                        d.append(std::to_string(dn));
                    }
                    if (pv) {
                        char nb[32];
                        for (size_t l = 0; l < L; ++l) { snprintf(nb, sizeof nb, "\t%.4f", pv[l]); d.append(nb); }
                    } else for (size_t l = 0; l < L; ++l) d.append("\tNA");
                    d.push_back('\n');
                }
                if (ents.size() == 1) {
                    out->kind.push_back(1);
                    out->tgt.push_back(ents[0].t);
                    continue;
                }
                if (cfg.PT && !cfg.PT->idx.empty()) {
                    rbuf.assign(recs[0].qname);
                    auto it = cfg.PT->idx.find(rbuf);
                    if (it != cfg.PT->idx.end()) pv = &cfg.PT->val[(size_t)it->second * L];
                }
                std::sort(ents.begin(), ents.end(), [](const Ent& x, const Ent& y) { return x.t < y.t; });
                const size_t n = ents.size();
                bool has_dn = false;
                tot.resize(n);
                if (NMW)
                    for (size_t i = 0; i < n; ++i) {
                        tot[i] = (ents[i].e0 < NM_INF ? ents[i].e0 : 0) + (ents[i].e1 < NM_INF ? ents[i].e1 : 0) - best;
                        if (tot[i]) has_dn = true;
                    }
                key.clear();
                key.append((const char*)&n, sizeof n);
                for (auto& e : ents) key.append((const char*)&e.t, sizeof e.t);
                if (!pv && !has_dn) {
                    out->kind.push_back(2);
                } else {
                    key.push_back(pv ? 1 : 0);
                    if (pv) {
                        for (auto& e : ents) {
                            double fct;
                            const auto& mem = M.POOLG[e.t];
                            if (mem.empty()) fct = M.tlin[e.t] >= 0 ? pv[M.tlin[e.t]] : 1.0 / (double)L;
                            else {
                                mbuf.clear();
                                for (uint32_t m : mem) mbuf.push_back(M.tlin[m] >= 0 ? pv[M.tlin[m]] : 1.0 / (double)L);
                                fct = np_sum(mbuf) / (double)mbuf.size();
                            }
                            fct = py_round4(fct);
                            key.append((const char*)&fct, sizeof fct);
                        }
                    }
                    key.push_back(has_dn ? 1 : 0);
                    if (has_dn) for (size_t i = 0; i < n; ++i) key.append((const char*)&tot[i], sizeof(int));
                    out->kind.push_back(3);
                }
                out->koff.push_back(out->keys.size());
                out->keys.append(key);
            }
            out->koff.push_back(out->keys.size());
            {
                std::lock_guard<std::mutex> lk(mu);
                results[b->id] = std::move(out);
            }
            cv_res.notify_all();
        }
    };
    std::vector<std::thread> ws;
    for (int i = 0; i < W; ++i) ws.emplace_back(worker, i);

    // ---------------------------------------------------------------- ordered merge
    std::string k;
    for (;;) {
        std::unique_ptr<BatchOut> r;
        {
            std::unique_lock<std::mutex> lk(mu);
            cv_res.wait(lk, [&]() { return results.count(merged) || (reader_done && merged == produced); });
            auto it = results.find(merged);
            if (it == results.end()) break;
            r = std::move(it->second);
            results.erase(it);
        }
        size_t ti = 0, ki = 0;
        for (uint8_t kd : r->kind) {
            if (kd == 1) { ev.single[r->tgt[ti++]] += 1.0; continue; }
            k.assign(r->keys, r->koff[ki], r->koff[ki + 1] - r->koff[ki]);
            ++ki;
            if (kd == 2) ev.eq.add(k, 1);
            else ev.eqx.add(k, 1);
        }
        ev.n_rec += r->n_rec; ev.n_frag += r->n_frag; ev.n_unal += r->n_unal;
        ev.nm_sum += r->nm_sum; ev.len_sum += r->len_sum;
        for (int32_t v : r->fl) { if (ev.fl_vals.size() >= 1000000) break; ev.fl_vals.push_back(v); }
        if (!r->gcb.empty()) {
            if (ev.gc_fg.empty()) ev.gc_fg.assign(101, 0.0);
            for (int i = 0; i <= 100; ++i) ev.gc_fg[i] += (double)r->gcb[i];
        }
        if (cfg.dbg_kmer && !r->dk.empty()) fwrite(r->dk.data(), 1, r->dk.size(), cfg.dbg_kmer);
        if (cfg.dbg_snp && !r->ds.empty()) fwrite(r->ds.data(), 1, r->ds.size(), cfg.dbg_snp);
        if (cfg.pvals_out && !r->pr.empty()) fwrite(r->pr.data(), 1, r->pr.size(), cfg.pvals_out);
        if (cfg.dbg_frag && !r->df.empty()) fwrite(r->df.data(), 1, r->df.size(), cfg.dbg_frag);
        {
            std::lock_guard<std::mutex> lk(mu);
            ++merged;
        }
        cv_space.notify_all();
    }
    reader.join();
    for (auto& t : ws) t.join();
    for (auto& a : accs) {
        for (size_t i = 0; i < a.im.size(); ++i) { ev.idx_mapped[i] += a.im[i]; ev.idx_unmapped[i] += a.iu[i]; }
        ev.idx_unmapped_star += a.star;
        ev.ls.add(a.ls);
    }
    if (cfg.LE) ev.n_pv_reads = ev.ls.pv_used;
}

}  // namespace lq
