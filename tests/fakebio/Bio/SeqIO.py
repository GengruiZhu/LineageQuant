"""Test stand-in for Bio.SeqIO.parse(path, "fasta") with Biopython's FASTA rules (SimpleFastaParser):
id = first word of the title; sequence = lines right-stripped, joined, spaces and \\r removed."""


class _Rec:
    def __init__(self, title, seq):
        self.id = title.split(None, 1)[0] if title else ""
        self.seq = seq


def parse(handle, fmt):
    assert fmt == "fasta"
    f = open(handle) if isinstance(handle, str) else handle
    title, lines = None, []
    for line in f:
        if line.startswith(">"):
            if title is not None:
                yield _Rec(title, "".join(lines).replace(" ", "").replace("\r", ""))
            title, lines = line[1:].rstrip(), []
        elif title is not None:
            lines.append(line.rstrip())
    if title is not None:
        yield _Rec(title, "".join(lines).replace(" ", "").replace("\r", ""))
