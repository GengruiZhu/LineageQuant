"""Test stand-in for cyvcf2 (plain / gzip VCF): VCF(path).samples, iteration over records with CHROM, POS, REF,
ALT and format(FIELD) -> int32 array (n_samples, max values) like htslib: missing = INT32_MIN, padding =
INT32_MIN + 1; None when the record's FORMAT lacks the field."""
import gzip
import numpy as np

MISSING, VEND = -2147483648, -2147483647


class _Rec:
    def __init__(self, f):
        self.CHROM, self.POS, self.REF = f[0], int(f[1]), f[3]
        self.ALT = [] if f[4] == "." else f[4].split(",")
        self._fmt = f[8].split(":") if len(f) > 8 else []
        self._s = f[9:]

    def format(self, field):
        if field not in self._fmt:
            return None
        j = self._fmt.index(field)
        vals = []
        for s in self._s:
            x = s.split(":")
            v = x[j] if j < len(x) else "."
            if v in (".", ""):
                vals.append([MISSING])
            else:
                vals.append([MISSING if t in (".", "") else int(t) for t in v.split(",")])
        m = max(len(v) for v in vals) if vals else 1
        return np.array([v + [VEND] * (m - len(v)) for v in vals], dtype=np.int32)


class VCF:
    def __init__(self, path):
        op = gzip.open if open(path, "rb").read(2) == b"\x1f\x8b" else open
        self._f = op(path, "rt")
        self.samples = []
        for line in self._f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                self.samples = line.rstrip("\n").split("\t")[9:]
                break

    def __iter__(self):
        for line in self._f:
            if line.strip():
                yield _Rec(line.rstrip("\n").split("\t"))
