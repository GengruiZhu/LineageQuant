"""Minimal stand-in for pysam, enough to run snipe_kmer.py / snipe_snp.py on a SAM file in tests.

AlignmentFile(path) accepts a SAM text file (name-grouped or not); fetch(ref) returns that reference's
records sorted by position, like an indexed coordinate-sorted BAM.
"""
_CACHE = {}


def _cig(s):
    out, n = [], ""
    for c in s:
        if c.isdigit():
            n += c
        else:
            out.append(("MIDNSHP=X".index(c), int(n)))
            n = ""
    return out


class AlignedSegment:
    __slots__ = ("query_name", "flag", "reference_name", "reference_start", "cigartuples", "_seq", "_tags")

    def __init__(self, q):
        self.query_name = q[0]
        self.flag = int(q[1])
        self.reference_name = q[2]
        self.reference_start = int(q[3]) - 1
        self.cigartuples = None if q[5] == "*" else _cig(q[5])
        self._seq = None if q[9] == "*" else q[9]
        self._tags = q[11:]

    @property
    def is_unmapped(self): return bool(self.flag & 4)

    @property
    def is_read1(self): return bool(self.flag & 64)

    @property
    def is_secondary(self): return bool(self.flag & 256)

    @property
    def query_sequence(self): return self._seq

    def get_tag(self, t):
        for x in self._tags:
            if x.startswith(t + ":"):
                v = x.split(":", 2)[2]
                return int(v) if x[len(t) + 1] == "i" else v
        raise KeyError(t)

    @property
    def reference_end(self):
        if self.cigartuples is None: return None
        return self.reference_start + sum(l for o, l in self.cigartuples if o in (0, 2, 3, 7, 8))

    def get_reference_positions(self, full_length=False):
        out, r = [], self.reference_start
        for o, l in self.cigartuples or []:
            if o in (0, 7, 8):
                out.extend(range(r, r + l)); r += l
            elif o in (1, 4):
                if full_length: out.extend([None] * l)
            elif o in (2, 3):
                r += l
        return out


class AlignmentFile:
    def __init__(self, path, mode="rb", check_sq=True):
        if path not in _CACHE:
            refs, by = [], {}
            with open(path) as f:
                for l in f:
                    if l.startswith("@SQ"):
                        refs.append(dict(x.split(":", 1) for x in l.rstrip("\n").split("\t")[1:])["SN"])
                    elif not l.startswith("@") and l.strip():
                        q = l.rstrip("\n").split("\t")
                        if q[2] != "*" and not int(q[1]) & 4:
                            by.setdefault(q[2], []).append((int(q[3]), len(by.get(q[2], [])), q))
            for k in by:
                by[k] = [AlignedSegment(x[2]) for x in sorted(by[k], key=lambda x: (x[0], x[1]))]
            _CACHE[path] = (tuple(refs), by)
        self.references, self._by = _CACHE[path]

    def fetch(self, tid=None, until_eof=False):
        return iter(self._by.get(tid, []))

    def close(self): pass
