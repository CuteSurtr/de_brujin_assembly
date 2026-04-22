"""Simulate short reads from a reference genome for pipeline validation."""

from __future__ import annotations

from typing import List, Optional

import numpy as np


def simulate_reads(
    reference: str,
    read_length: int = 100,
    coverage: float = 30.0,
    error_rate: float = 0.0,
    rng: Optional[np.random.Generator] = None,
) -> List[str]:
    """Generate short reads covering the reference at the target coverage.

    Reads are sampled uniformly at random from the reference. Errors, if
    any, are introduced as random substitutions.
    """
    if rng is None:
        rng = np.random.default_rng()
    n = len(reference)
    n_reads = int(round(coverage * n / read_length))
    reads: List[str] = []
    bases = "ACGT"
    for _ in range(n_reads):
        start = int(rng.integers(0, max(n - read_length + 1, 1)))
        read = reference[start : start + read_length]
        if error_rate > 0 and read:
            r_arr = bytearray(read, "ascii")
            for i in range(len(r_arr)):
                if rng.random() < error_rate:
                    other = bases.replace(chr(r_arr[i]), "")
                    r_arr[i] = ord(other[int(rng.integers(0, 3))])
            read = r_arr.decode("ascii")
        reads.append(read)
    return reads
