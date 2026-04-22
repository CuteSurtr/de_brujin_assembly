from __future__ import annotations

from pathlib import Path
from typing import Iterator, List, Tuple

DNA_ALPHABET = "ACGT"
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def canonical(kmer: str) -> str:
    rc = revcomp(kmer)
    return kmer if kmer <= rc else rc


def read_fasta(path) -> List[Tuple[str, str]]:
    text = Path(path).read_text()
    out: List[Tuple[str, str]] = []
    header = None
    parts: List[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                out.append((header, "".join(parts).upper()))
            header = line[1:].strip()
            parts = []
        else:
            parts.append(line.strip())
    if header is not None:
        out.append((header, "".join(parts).upper()))
    return out


def read_fastq(path) -> List[Tuple[str, str]]:
    text = Path(path).read_text()
    lines = text.splitlines()
    out: List[Tuple[str, str]] = []
    for i in range(0, len(lines), 4):
        if i + 1 >= len(lines):
            break
        name = lines[i].lstrip("@").strip()
        seq = lines[i + 1].strip().upper()
        out.append((name, seq))
    return out


def write_fasta(path, records: List[Tuple[str, str]]) -> None:
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def kmers(seq: str, k: int) -> Iterator[str]:
    for i in range(len(seq) - k + 1):
        km = seq[i : i + k]
        if "N" not in km:
            yield km
