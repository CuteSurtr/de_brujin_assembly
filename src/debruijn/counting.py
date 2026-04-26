from __future__ import annotations

from collections import Counter
from typing import Dict, Iterable, List

import numpy as np

from .utils import canonical, kmers, revcomp


def count_kmers(reads: Iterable[str], k: int, canonicalize: bool = True) -> Counter:
    """Exact k-mer count over a set of reads. If ``canonicalize`` is True,
    each k-mer and its reverse complement are counted together under the
    canonical (lexicographically smaller) form, which is standard for
    double stranded DNA."""
    c: Counter = Counter()
    for read in reads:
        for km in kmers(read, k):
            if canonicalize:
                c[canonical(km)] += 1
            else:
                c[km] += 1
    return c


def kmer_histogram(counts: Counter, max_count: int = 100) -> np.ndarray:
    """Histogram of how many distinct k-mers occur exactly 1, 2, ..., max_count times.
    Useful for choosing a coverage cutoff -- the first valley between the error peak
    near 1 and the genomic peak near mean_coverage is the natural threshold."""
    h = np.zeros(max_count + 1, dtype=int)
    for cnt in counts.values():
        h[min(cnt, max_count)] += 1
    return h


class CountMinSketch:
    """Approximate k-mer counter using a small count-min sketch. Gives a
    probabilistic over estimate of any given k-mer's count in sub-linear
    space -- useful when the distinct k-mer set is too large for a dict."""

    def __init__(self, width: int = 1 << 20, depth: int = 4, seed: int = 0):
        self.width = width
        self.depth = depth
        self.table = np.zeros((depth, width), dtype=np.int32)
        rng = np.random.default_rng(seed)
        self._seeds = rng.integers(1, 2 ** 31 - 1, size=depth, dtype=np.int64)

    def _indices(self, item: str):
        # Simple per row hash: (hash(item) ^ seed_i) mod width
        h = hash(item)
        return [(h ^ int(s)) % self.width for s in self._seeds]

    def add(self, item: str, count: int = 1) -> None:
        for i, idx in enumerate(self._indices(item)):
            self.table[i, idx] += count

    def query(self, item: str) -> int:
        return int(min(self.table[i, idx] for i, idx in enumerate(self._indices(item))))
