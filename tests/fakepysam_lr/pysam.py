"""Stand-in for pysam for the long-read tests (extract_prior.py / calibrate_snp.py).

AlignmentFile(path) reads a SAM text file; fetch(until_eof=True) yields records in file order
(including unmapped ones), like iterating an unindexed BAM.
"""


def _cig(s):
    out, n = [], ""
    for c in s:
        if c.isdigit():
            n += c
        else:
            out.append(("MIDNSHP=X".index(c), int(n)))
            n = ""
    return out


class _Header:
    def __init__(self, hd, sq):
        self._hd, self._sq = hd, sq

    def to_dict(self):
        return {"HD": dict(self._hd), "SQ": [{"SN": s} for s in self._sq]}


class AlignedSegment:
    def __init__(self, q, r2i):
        self.query_name = q[0]
        self.flag = int(q[1])
        self.reference_name = q[2]
        self.reference_id = r2i.get(q[2], -1)
        self.reference_start = int(q[3]) - 1
        self.cigartuples = None if q[5] == "*" else _cig(q[5])
        self._seq = None if q[9] == "*" else q[9]
        self._tags = q[11:]

    is_unmapped = property(lambda s: bool(s.flag & 4))
    is_reverse = property(lambda s: bool(s.flag & 16))
    is_secondary = property(lambda s: bool(s.flag & 256))
    is_qcfail = property(lambda s: bool(s.flag & 512))
    is_duplicate = property(lambda s: bool(s.flag & 1024))
    is_supplementary = property(lambda s: bool(s.flag & 2048))
    query_sequence = property(lambda s: s._seq)

    def has_tag(self, t):
        return any(x.startswith(t + ":") for x in self._tags)

    def get_tag(self, t):
        for x in self._tags:
            if x.startswith(t + ":"):
                v = x.split(":", 2)[2]
                return int(v) if x[len(t) + 1] == "i" else v
        raise KeyError(t)

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
        hd, sq, self._recs = {}, [], []
        with open(path) as f:
            lines = f.read().splitlines()
        for l in lines:
            if l.startswith("@HD"):
                hd = dict(x.split(":", 1) for x in l.split("\t")[1:])
            elif l.startswith("@SQ"):
                sq.append(dict(x.split(":", 1) for x in l.split("\t")[1:])["SN"])
        self.references = tuple(sq)
        self._r2i = {s: i for i, s in enumerate(sq)}
        self._lines = [l.split("\t") for l in lines if l and not l.startswith("@")]
        self.header = _Header(hd, sq)

    def get_reference_name(self, i):
        return self.references[i]

    def fetch(self, *a, until_eof=False, **k):
        assert until_eof, "test BAMs have no index"
        return (AlignedSegment(q, self._r2i) for q in self._lines)

    def close(self): pass

    def __enter__(self): return self

    def __exit__(self, *a): self.close()
