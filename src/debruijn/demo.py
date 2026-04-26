"""End to end demonstration of de Bruijn graph assembly on lambda phage.

The demo simulates short reads from the lambda phage reference genome
(48502 bp, NCBI NC_001416.1), builds a de Bruijn graph, cleans it via
tip removal and low coverage edge filtering, extracts contigs as
maximal non branching paths, and reports assembly statistics. It
sweeps coverage depth and k to illustrate the classical k vs coverage
tradeoff documented in Zerbino and Birney 2008, Figure 3.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

import numpy as np

from . import (
    DeBruijnGraph,
    assembly_stats,
    count_kmers,
    extract_contigs,
    find_eulerian_path,
    kmer_histogram,
    n50,
    path_to_sequence,
    read_fasta,
    simulate_reads,
)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False


HERE = Path(__file__).resolve().parents[2]
DATA = HERE / "data"
RESULTS = HERE / "results"


def section(title: str) -> None:
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74)


def assemble(
    reads: List[str],
    k: int,
    min_coverage: int = 2,
    tip_multiple: int = 2,
) -> Tuple[List[str], dict]:
    """Build, clean, and extract contigs from a read set."""
    g = DeBruijnGraph.from_reads(reads, k=k)
    g.remove_low_coverage_edges(min_coverage=min_coverage)
    g.remove_tips(max_tip_length=tip_multiple * k)
    contigs = extract_contigs(g)
    return contigs, assembly_stats(contigs)


def main() -> None:
    RESULTS.mkdir(exist_ok=True)

    # ---- load reference ---------------------------------------------------
    fasta_path = DATA / "lambda_phage.fasta"
    if not fasta_path.exists():
        print(f"ERROR: could not find {fasta_path}")
        print("Download with:")
        print('  curl -sL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
              'efetch.fcgi?db=nuccore&id=NC_001416.1&rettype=fasta" '
              '-o data/lambda_phage.fasta')
        return
    records = read_fasta(fasta_path)
    name, reference = records[0]
    reference = reference.upper().replace("N", "")
    print(f"loaded {name.split()[0]}: {len(reference):,} bp")

    # ---- 1. simulate error free reads at modest coverage ------------------
    section("1. Error free simulation")
    rng = np.random.default_rng(42)
    reads = simulate_reads(reference, read_length=100, coverage=20.0,
                           error_rate=0.0, rng=rng)
    print(f"  simulated {len(reads):,} reads of length 100")
    print(f"  effective coverage = {sum(len(r) for r in reads) / len(reference):.1f}x")

    contigs, stats = assemble(reads, k=31, min_coverage=2)
    print(f"  k=31, min_cov=2 -> {stats['n_contigs']:,} contigs, "
          f"total {stats['total_length']:,} bp, "
          f"N50 = {stats['n50']:,} bp, longest = {stats['longest']:,} bp")
    print(f"  recovered = {stats['total_length'] / len(reference) * 100:.1f}% of reference")

    # ---- 2. sweep k at fixed coverage -------------------------------------
    section("2. Choice of k at fixed coverage (30x, error free)")
    rng = np.random.default_rng(0)
    reads = simulate_reads(reference, read_length=100, coverage=30.0,
                           error_rate=0.0, rng=rng)
    print(f"  {'k':>3}  {'n_contigs':>10}  {'N50':>10}  {'longest':>10}  {'recovered':>10}")
    k_sweep_results = []
    for k in [15, 21, 25, 31, 41, 51, 71]:
        if k >= 100:
            continue
        _, stats = assemble(reads, k=k, min_coverage=2)
        rec = stats['total_length'] / len(reference)
        k_sweep_results.append((k, stats))
        print(f"  {k:>3}  {stats['n_contigs']:>10,}  {stats['n50']:>10,}  "
              f"{stats['longest']:>10,}  {rec * 100:>9.1f}%")

    # ---- 3. sweep coverage at fixed k -------------------------------------
    section("3. Coverage sweep at k=31 (error free)")
    print(f"  {'cov':>5}  {'n_contigs':>10}  {'N50':>10}  {'longest':>10}")
    cov_sweep_results = []
    for cov in [5, 10, 20, 30, 50]:
        rng = np.random.default_rng(cov)
        reads = simulate_reads(reference, read_length=100, coverage=float(cov),
                               error_rate=0.0, rng=rng)
        _, stats = assemble(reads, k=31, min_coverage=max(2, cov // 10))
        cov_sweep_results.append((cov, stats))
        print(f"  {cov:>4}x  {stats['n_contigs']:>10,}  {stats['n50']:>10,}  "
              f"{stats['longest']:>10,}")

    # ---- 4. effect of sequencing errors -----------------------------------
    section("4. Effect of sequencing errors at k=31, 30x coverage")
    print(f"  {'error':>7}  {'n_contigs':>10}  {'N50':>10}  {'longest':>10}")
    err_sweep_results = []
    for err in [0.0, 0.001, 0.005, 0.01, 0.02]:
        rng = np.random.default_rng(int(err * 10000) + 1)
        reads = simulate_reads(reference, read_length=100, coverage=30.0,
                               error_rate=err, rng=rng)
        # With errors, raise the coverage cutoff to filter spurious k-mers.
        min_cov = 3 if err == 0 else max(3, int(30 * err * 10))
        _, stats = assemble(reads, k=31, min_coverage=min_cov)
        err_sweep_results.append((err, stats))
        print(f"  {err * 100:>6.2f}%  {stats['n_contigs']:>10,}  "
              f"{stats['n50']:>10,}  {stats['longest']:>10,}")

    # ---- 5. k-mer spectrum inspection -------------------------------------
    section("5. k-mer spectrum at k=21, 30x error free coverage")
    rng = np.random.default_rng(99)
    reads = simulate_reads(reference, read_length=100, coverage=30.0,
                           error_rate=0.005, rng=rng)
    counts = count_kmers(reads, k=21, canonicalize=True)
    hist = kmer_histogram(counts, max_count=80)
    total_distinct = int(hist.sum())
    print(f"  {total_distinct:,} distinct 21-mers observed")
    print("  coverage histogram (first 20 bins):")
    h_max = int(hist.max())
    for c in range(1, 21):
        bar = "#" * int(int(hist[c]) / max(h_max, 1) * 40)
        print(f"    cov={c:>3}: {int(hist[c]):>8,}  {bar}")
    print("  (the peak near expected coverage ~ 30 represents genuine k-mers;")
    print("   the tall spike at cov=1 is sequencing error k-mers.)")

    # ---- 6. Eulerian reconstruction on small substring --------------------
    section("6. Full Eulerian reconstruction on a 2 kb substring (no errors)")
    sub = reference[10000:12000]
    rng = np.random.default_rng(7)
    reads_sub = simulate_reads(sub, read_length=80, coverage=25.0,
                               error_rate=0.0, rng=rng)
    g_sub = DeBruijnGraph.from_reads(reads_sub, k=25)
    print(f"  graph: {len(g_sub.nodes()):,} nodes, "
          f"{g_sub.distinct_edge_count():,} distinct edges "
          f"({g_sub.edge_count():,} with multiplicity)")
    path = find_eulerian_path(g_sub)
    if path is not None:
        recon = path_to_sequence(path)
        from .utils import kmers as _kmers
        ref_kmers = set(_kmers(sub, 25))
        recon_kmers = set(_kmers(recon, 25))
        shared = ref_kmers & recon_kmers
        print(f"  Eulerian path found: length {len(path):,} nodes")
        print(f"  reference k-mers: {len(ref_kmers):,}, "
              f"reconstructed k-mers: {len(recon_kmers):,}, "
              f"shared: {len(shared):,}")
        print(f"  reference k-mer recovery: {len(shared) / len(ref_kmers) * 100:.1f}%")
    else:
        print("  (no Eulerian path -- fell back to contig extraction)")

    # ---- 7. figures -------------------------------------------------------
    if HAVE_MPL:
        section("7. Figures -> results/")
        fig, axs = plt.subplots(2, 2, figsize=(13, 9))

        # (a) N50 vs k
        ks = [k for k, _ in k_sweep_results]
        n50s = [s["n50"] for _, s in k_sweep_results]
        longests = [s["longest"] for _, s in k_sweep_results]
        axs[0][0].plot(ks, n50s, "o-", label="N50", color="#1f77b4")
        axs[0][0].plot(ks, longests, "s--", label="longest contig", color="#ff7f0e")
        axs[0][0].set_xlabel("k")
        axs[0][0].set_ylabel("bp")
        axs[0][0].set_title("Contig quality vs. k (30x error-free)")
        axs[0][0].legend()
        axs[0][0].grid(True, alpha=0.3)

        # (b) N50 vs coverage
        covs = [c for c, _ in cov_sweep_results]
        c_n50s = [s["n50"] for _, s in cov_sweep_results]
        c_long = [s["longest"] for _, s in cov_sweep_results]
        axs[0][1].plot(covs, c_n50s, "o-", label="N50", color="#2ca02c")
        axs[0][1].plot(covs, c_long, "s--", label="longest", color="#d62728")
        axs[0][1].set_xlabel("coverage (x)")
        axs[0][1].set_ylabel("bp")
        axs[0][1].set_title("Contig quality vs. coverage (k=31)")
        axs[0][1].legend()
        axs[0][1].grid(True, alpha=0.3)

        # (c) N contigs vs error rate
        errs = [e * 100 for e, _ in err_sweep_results]
        ns = [s["n_contigs"] for _, s in err_sweep_results]
        axs[1][0].plot(errs, ns, "o-", color="#9467bd")
        axs[1][0].set_xlabel("sequencing error rate (%)")
        axs[1][0].set_ylabel("number of contigs")
        axs[1][0].set_title("Fragmentation vs. error rate (k=31, 30x)")
        axs[1][0].grid(True, alpha=0.3)

        # (d) k-mer spectrum
        bins = list(range(1, len(hist)))
        heights = [int(hist[b]) for b in bins]
        axs[1][1].bar(bins, heights, color="#8c564b", edgecolor="black", linewidth=0.3)
        axs[1][1].set_yscale("log")
        axs[1][1].set_xlabel("k-mer coverage")
        axs[1][1].set_ylabel("number of distinct k-mers (log scale)")
        axs[1][1].set_title("k-mer spectrum (k=21, 30x with 0.5% error)")
        axs[1][1].grid(True, alpha=0.3, which="both")

        fig.tight_layout()
        fig.savefig(RESULTS / "debruijn_demo.png", dpi=140, bbox_inches="tight")
        print(f"  wrote {RESULTS / 'debruijn_demo.png'}")


if __name__ == "__main__":
    main()
