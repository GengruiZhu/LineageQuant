"""Test stand-in for Bio.AlignIO.read(path, "fasta"): records of equal length, else ValueError."""
from Bio import SeqIO


class _Aln(list):
    def get_alignment_length(self):
        return len(self[0].seq)


def read(handle, fmt):
    recs = list(SeqIO.parse(handle, fmt))
    if not recs:
        raise ValueError("No records found in handle")
    if len(set(len(r.seq) for r in recs)) != 1:
        raise ValueError("Sequences must all be the same length")
    return _Aln(recs)
