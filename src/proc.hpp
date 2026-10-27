// LineageQuant 1.1.0 -- running external tools (bowtie2-build, salmon, kmc, kmc_tools, gffread, mafft, samtools)
#pragma once
#include "common.hpp"

#include <csignal>
#include <fcntl.h>
#include <poll.h>
#include <sys/wait.h>

namespace lq {

// PATH lookup like shutil.which(); "" if not found
inline std::string which(const std::string& tool) {
    auto exe = [](const std::string& p) {
        struct stat st;
        return stat(p.c_str(), &st) == 0 && S_ISREG(st.st_mode) && access(p.c_str(), X_OK) == 0;
    };
    if (tool.find('/') != std::string::npos) return exe(tool) ? tool : std::string();
    const char* path = getenv("PATH");
    std::string P = path ? path : "/usr/bin:/bin";
    size_t b = 0;
    for (;;) {
        size_t e = P.find(':', b);
        std::string d = P.substr(b, e == std::string::npos ? std::string::npos : e - b);
        std::string c = (d.empty() ? std::string(".") : d) + "/" + tool;
        if (exe(c)) return c;
        if (e == std::string::npos) return {};
        b = e + 1;
    }
}

inline std::string shell_quote(const std::string& s) {
    bool plain = !s.empty();
    for (char c : s)
        if (!(isalnum((unsigned char)c) || strchr("@%+=:,./-_", c))) { plain = false; break; }
    if (plain) return s;
    std::string o = "'";
    for (char c : s) o += c == '\'' ? std::string("'\"'\"'") : std::string(1, c);
    return o + "'";
}
inline std::string join_cmd(const std::vector<std::string>& v) {
    std::string s;
    for (auto& x : v) s += (s.empty() ? "" : " ") + shell_quote(x);
    return s;
}

// Run argv[0] (looked up on PATH) without a shell.
//   out / err: "" = inherit, "-" = /dev/null, "|" = capture into *cap (stdout only), else a file path
//   (err may be "&1": same as stdout); timeout_s > 0 kills the process after that many seconds.
// Returns the exit status (0 = ok), 124 on timeout, 127 if the program cannot be started.
struct RunResult {
    int status = 0;
    bool timed_out = false;
};
inline RunResult run_proc(const std::vector<std::string>& argv, const std::string& out = "", const std::string& err = "",
                          std::string* cap = nullptr, int timeout_s = 0) {
    RunResult R;
    std::string prog = which(argv[0]);
    if (prog.empty()) { R.status = 127; return R; }
    std::vector<char*> av;
    for (auto& s : argv) av.push_back(const_cast<char*>(s.c_str()));
    av.push_back(nullptr);
    int pfd[2] = {-1, -1};
    if (out == "|" && pipe(pfd) != 0) die("pipe() failed");
    fflush(stdout);
    fflush(stderr);
    pid_t pid = fork();
    if (pid < 0) die("fork() failed");
    if (pid == 0) {
        setpgid(0, 0);
        signal(SIGPIPE, SIG_DFL);
        auto redirect = [](const std::string& w, int fd, int cap_fd) {
            if (w.empty()) return;
            int t;
            if (w == "-") t = open("/dev/null", O_WRONLY);
            else if (w == "|") t = cap_fd;
            else if (w == "&1") t = dup(1);
            else t = open(w.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
            if (t < 0) _exit(126);
            dup2(t, fd);
            if (t != fd) close(t);
        };
        if (pfd[0] >= 0) close(pfd[0]);
        redirect(out, 1, pfd[1]);
        redirect(err, 2, -1);
        execv(prog.c_str(), av.data());
        _exit(127);
    }
    setpgid(pid, pid);
    if (pfd[1] >= 0) close(pfd[1]);
    const double t0 = wall();
    if (pfd[0] >= 0) {
        char buf[1 << 16];
        for (;;) {
            int wait_ms = -1;
            if (timeout_s > 0) {
                double left = timeout_s - (wall() - t0);
                if (left <= 0) { R.timed_out = true; break; }
                wait_ms = (int)(left * 1000) + 1;
            }
            struct pollfd q = {pfd[0], POLLIN, 0};
            int r = poll(&q, 1, wait_ms);
            if (r < 0) { if (errno == EINTR) continue; break; }
            if (r == 0) { R.timed_out = true; break; }
            ssize_t n = read(pfd[0], buf, sizeof buf);
            if (n < 0) { if (errno == EINTR) continue; break; }
            if (n == 0) break;
            if (cap) cap->append(buf, (size_t)n);
        }
        close(pfd[0]);
    }
    int st = 0;
    for (;;) {
        pid_t w = waitpid(pid, &st, timeout_s > 0 ? WNOHANG : 0);
        if (w == pid) break;
        if (w < 0 && errno != EINTR) { st = 0; break; }
        if (timeout_s > 0) {
            if (R.timed_out || wall() - t0 > timeout_s) {
                R.timed_out = true;
                kill(-pid, SIGKILL);
                waitpid(pid, &st, 0);
                break;
            }
            usleep(20000);
        }
    }
    if (R.timed_out) { R.status = 124; return R; }
    R.status = WIFEXITED(st) ? WEXITSTATUS(st) : 128 + (WIFSIGNALED(st) ? WTERMSIG(st) : 0);
    return R;
}

// run, die on failure (subprocess.run(check=True))
inline void run_checked(const std::vector<std::string>& argv, const std::string& out = "", const std::string& err = "") {
    RunResult r = run_proc(argv, out, err);
    if (r.status == 127) die("cannot run '" + argv[0] + "' (not on PATH?)");
    if (r.status != 0) die("command failed (exit " + std::to_string(r.status) + "): " + join_cmd(argv));
}

inline long long file_size(const std::string& p) {
    struct stat st;
    return stat(p.c_str(), &st) == 0 ? (long long)st.st_size : -1;
}
inline bool nonempty(const std::string& p) { return file_size(p) > 0; }
inline void touch_done(const std::string& p) {
    FILE* f = fopen(p.c_str(), "w");
    if (!f) die("cannot write " + p);
    fputs("DONE", f);
    fclose(f);
}
}  // namespace lq
