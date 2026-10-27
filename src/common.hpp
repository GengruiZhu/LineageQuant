// lqcore -- LineageQuant core in C++
// common.hpp: logging, formatting, line reading, numpy-compatible arithmetic helpers
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <string_view>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <zlib.h>

namespace lq {

// ------------------------------------------------------------------ logging
inline std::string now_str() {
    time_t t = time(nullptr);
    struct tm tm;
    localtime_r(&t, &tm);
    char b[32];
    strftime(b, sizeof b, "%Y-%m-%d %H:%M:%S", &tm);
    return b;
}
inline void logm(const std::string& m) {
    std::printf("[%s] %s\n", now_str().c_str(), m.c_str());
    std::fflush(stdout);
}
[[noreturn]] inline void die(const std::string& m) {
    std::fflush(stdout);
    std::fprintf(stderr, "[%s] [FATAL] %s\n", now_str().c_str(), m.c_str());
    std::exit(1);
}
inline double wall() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// ------------------------------------------------------------------ formatting
inline std::string sfmt(const char* f, ...) __attribute__((format(printf, 1, 2)));
inline std::string sfmt(const char* f, ...) {
    va_list ap;
    va_start(ap, f);
    char buf[1024];
    int n = vsnprintf(buf, sizeof buf, f, ap);
    va_end(ap);
    if (n < 0) return {};
    if (n < (int)sizeof buf) return std::string(buf, n);
    std::string s(n, '\0');
    va_start(ap, f);
    vsnprintf(&s[0], n + 1, f, ap);
    va_end(ap);
    return s;
}
// Python f"{v:,}" for integers
inline std::string commas(long long v) {
    bool neg = v < 0;
    unsigned long long u = neg ? (unsigned long long)(-v) : (unsigned long long)v;
    std::string s = std::to_string(u), o;
    int c = 0;
    for (int i = (int)s.size() - 1; i >= 0; --i) {
        o.push_back(s[i]);
        if (++c % 3 == 0 && i > 0) o.push_back(',');
    }
    if (neg) o.push_back('-');
    std::reverse(o.begin(), o.end());
    return o;
}
// Python f"{x:,.0f}"
inline std::string commas_f0(double x) {
    std::string s = sfmt("%.0f", x);
    bool neg = !s.empty() && s[0] == '-';
    std::string d = neg ? s.substr(1) : s, o;
    int c = 0;
    for (int i = (int)d.size() - 1; i >= 0; --i) {
        o.push_back(d[i]);
        if (++c % 3 == 0 && i > 0) o.push_back(',');
    }
    if (neg) o.push_back('-');
    std::reverse(o.begin(), o.end());
    return o;
}
// Python repr-like short float for option echo (0.001 -> "0.001")
inline std::string pyfloat(double x) {
    for (int p = 1; p <= 17; ++p) {
        std::string s = sfmt("%.*g", p, x);
        if (strtod(s.c_str(), nullptr) == x) {
            if (s.find_first_of(".eEn") == std::string::npos) s += ".0";
            return s;
        }
    }
    return sfmt("%.17g", x);
}

// Python repr(float)
inline std::string py_repr(double x) {
    if (std::isnan(x)) return "nan";
    if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
    if (x == 0) return std::signbit(x) ? "-0.0" : "0.0";
    char buf[64];
    for (int p = 1; p <= 17; ++p) {
        snprintf(buf, sizeof buf, "%.*e", p - 1, x);
        if (strtod(buf, nullptr) == x) break;
    }
    std::string s(buf);
    bool neg = s[0] == '-';
    if (neg) s.erase(0, 1);
    size_t e = s.find('e');
    int ex = std::atoi(s.c_str() + e + 1);
    std::string d;
    for (size_t i = 0; i < e; ++i) if (isdigit((unsigned char)s[i])) d += s[i];
    while (d.size() > 1 && d.back() == '0') d.pop_back();
    const int n = (int)d.size();
    std::string out;
    if (ex >= -4 && ex < 16) {
        if (ex >= 0) out = n <= ex + 1 ? d + std::string(ex + 1 - n, '0') + ".0" : d.substr(0, ex + 1) + "." + d.substr(ex + 1);
        else out = "0." + std::string(-ex - 1, '0') + d;
    } else {
        out = d.substr(0, 1);
        if (n > 1) out += "." + d.substr(1);
        out += sfmt("e%c%02d", ex < 0 ? '-' : '+', std::abs(ex));
    }
    return neg ? "-" + out : out;
}

// ------------------------------------------------------------------ numpy-compatible arithmetic
// np.add.reduce on a contiguous float64 array (pairwise summation, PW_BLOCKSIZE = 128)
inline double np_pairwise(const double* a, size_t n) {
    if (n < 8) {
        double res = 0.;
        for (size_t i = 0; i < n; ++i) res += a[i];
        return res;
    }
    if (n <= 128) {
        double r[8];
        for (int j = 0; j < 8; ++j) r[j] = a[j];
        size_t i;
        for (i = 8; i < n - (n % 8); i += 8)
            for (int j = 0; j < 8; ++j) r[j] += a[i + j];
        double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; ++i) res += a[i];
        return res;
    }
    size_t n2 = n / 2;
    n2 -= n2 % 8;
    return np_pairwise(a, n2) + np_pairwise(a + n2, n - n2);
}
// np.add.reduce / ndarray.sum() of a float64 array: numpy reduces in buffer-sized blocks of 8192 elements
// (pairwise inside a block, left fold over blocks)
inline double np_add_reduce(const double* a, size_t n) {
    if (n <= 8192) return np_pairwise(a, n);
    double r = 0.0;
    for (size_t i = 0; i < n; i += 8192) r += np_pairwise(a + i, std::min<size_t>(8192, n - i));
    return r;
}
inline double np_sum(const std::vector<double>& v) { return np_add_reduce(v.data(), v.size()); }
// Python round(x, 4): correctly rounded decimal conversion, then back to double
inline double py_round(double x, int nd) {
    char b[400];
    snprintf(b, sizeof b, "%.*f", nd, x);
    return strtod(b, nullptr);
}
inline double py_round4(double x) {
    char b[64];
    snprintf(b, sizeof b, "%.4f", x);
    return strtod(b, nullptr);
}

// ------------------------------------------------------------------ strings
inline std::vector<std::string_view> split_tab(std::string_view s) {
    std::vector<std::string_view> out;
    size_t b = 0;
    for (;;) {
        size_t e = s.find('\t', b);
        if (e == std::string_view::npos) {
            out.emplace_back(s.substr(b));
            return out;
        }
        out.emplace_back(s.substr(b, e - b));
        b = e + 1;
    }
}
// Python str.split() (runs of whitespace)
inline std::vector<std::string_view> split_ws(std::string_view s) {
    std::vector<std::string_view> out;
    size_t i = 0, n = s.size();
    while (i < n) {
        while (i < n && isspace((unsigned char)s[i])) ++i;
        if (i >= n) break;
        size_t j = i;
        while (j < n && !isspace((unsigned char)s[j])) ++j;
        out.emplace_back(s.substr(i, j - i));
        i = j;
    }
    return out;
}
// Python s.rsplit(".", 1)[0]
inline std::string_view rsplit_dot(std::string_view s) {
    size_t p = s.rfind('.');
    return p == std::string_view::npos ? s : s.substr(0, p);
}
// re.sub(r'\.t?\d+$', '', s.split()[0])
inline std::string strip_iso(std::string_view g) {
    auto ws = split_ws(g);
    std::string_view s = ws.empty() ? std::string_view() : ws[0];
    size_t n = s.size(), i = n;
    while (i > 0 && isdigit((unsigned char)s[i - 1])) --i;
    if (i == n) return std::string(s);           // no trailing digits
    size_t j = i;
    if (j > 0 && s[j - 1] == 't') --j;
    if (j > 0 && s[j - 1] == '.') return std::string(s.substr(0, j - 1));
    if (i > 0 && s[i - 1] == '.') return std::string(s.substr(0, i - 1));   // ".\d+" with a preceding 't' that is not part of the suffix
    return std::string(s);
}
inline std::string upper(std::string_view s) {
    std::string o(s);
    for (auto& c : o) c = (char)toupper((unsigned char)c);
    return o;
}
// Python float(); returns false on failure (ValueError)
inline bool py_float(std::string_view s, double& v) {
    std::string t(s);
    size_t a = 0, b = t.size();
    while (a < b && isspace((unsigned char)t[a])) ++a;
    while (b > a && isspace((unsigned char)t[b - 1])) --b;
    if (a == b) return false;
    t = t.substr(a, b - a);
    char* end = nullptr;
    v = strtod(t.c_str(), &end);
    return end && *end == '\0';
}
inline bool file_exists(const std::string& p) { return access(p.c_str(), F_OK) == 0; }
// mkdir -p
inline void mkdirs(const std::string& d) {
    if (d.empty()) return;
    for (size_t i = 1; i <= d.size(); ++i)
        if (i == d.size() || d[i] == '/') mkdir(d.substr(0, i).c_str(), 0775);
}
// integer field of a SAM line (like int(); no allocation)
inline long long parse_ll(std::string_view s) {
    size_t i = 0;
    bool neg = false;
    if (i < s.size() && (s[i] == '-' || s[i] == '+')) neg = s[i++] == '-';
    long long v = 0;
    for (; i < s.size() && s[i] >= '0' && s[i] <= '9'; ++i) v = v * 10 + (s[i] - '0');
    return neg ? -v : v;
}

// ------------------------------------------------------------------ line reader (plain or gzip, file / stdin / fd)
class LineReader {
    gzFile gz_ = nullptr;
    std::vector<char> buf_;
    size_t beg_ = 0, end_ = 0;
    bool eof_ = false;

public:
    explicit LineReader(const std::string& path) {
        gz_ = (path == "-") ? gzdopen(dup(0), "rb") : gzopen(path.c_str(), "rb");
        if (!gz_) die("cannot open " + path);
        init();
    }
    explicit LineReader(int fd) {
        gz_ = gzdopen(fd, "rb");
        if (!gz_) die("cannot open file descriptor");
        init();
    }
    ~LineReader() {
        if (gz_) gzclose(gz_);
    }
    LineReader(const LineReader&) = delete;
    LineReader& operator=(const LineReader&) = delete;
    void init() {
        gzbuffer(gz_, 1 << 20);
        buf_.resize(1 << 24);
    }
    // line without the trailing "\n" (and "\r", like Python's universal newlines); valid until next call
    bool next(std::string_view& line) {
        for (;;) {
            char* base = buf_.data();
            char* p = (char*)memchr(base + beg_, '\n', end_ - beg_);
            if (p) {
                size_t len = p - (base + beg_);
                if (len && base[beg_ + len - 1] == '\r') --len;
                line = std::string_view(base + beg_, len);
                beg_ = (p - base) + 1;
                return true;
            }
            if (eof_) {
                if (beg_ < end_) {
                    size_t len = end_ - beg_;
                    if (len && base[beg_ + len - 1] == '\r') --len;
                    line = std::string_view(base + beg_, len);
                    beg_ = end_;
                    return true;
                }
                return false;
            }
            if (beg_ > 0) {
                memmove(base, base + beg_, end_ - beg_);
                end_ -= beg_;
                beg_ = 0;
            }
            if (end_ == buf_.size()) buf_.resize(buf_.size() * 2);
            size_t want = std::min<size_t>(buf_.size() - end_, (size_t)1 << 30);
            int n = gzread(gz_, buf_.data() + end_, (unsigned)want);
            if (n < 0) {
                int err;
                die(std::string("read error: ") + gzerror(gz_, &err));
            }
            if (n == 0) eof_ = true;
            else end_ += (size_t)n;
        }
    }
};

// ------------------------------------------------------------------ 128-bit sequence hash (identity of transcript sequences)
struct Hash128 {
    uint64_t a, b;
    bool operator==(const Hash128& o) const { return a == o.a && b == o.b; }
};
struct Hash128H {
    size_t operator()(const Hash128& h) const { return (size_t)(h.a ^ (h.b * 0x9E3779B97F4A7C15ULL)); }
};
class SeqHasher {
    uint64_t a_ = 0xcbf29ce484222325ULL, b_ = 0x84222325cbf29ce4ULL;

public:
    void add(std::string_view s) {
        for (unsigned char c : s) {
            c = (unsigned char)toupper(c);
            a_ = (a_ ^ c) * 0x100000001b3ULL;
            b_ = (b_ + c + 0x9E3779B97F4A7C15ULL) * 0xff51afd7ed558ccdULL;
            b_ ^= b_ >> 29;
        }
    }
    Hash128 get() const { return {a_, b_}; }
};

}  // namespace lq
