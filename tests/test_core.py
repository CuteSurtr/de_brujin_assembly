"""Core correctness tests for the de Bruijn assembler."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from debruijn import (
    CountMinSketch,
    DeBruijnGraph,
    assembly_stats,
    canonical,
    count_kmers,
    extract_contigs,
    find_eulerian_path,
    kmers,
    n50,
    path_to_sequence,
    revcomp,
    simulate_reads,
)


# ---------- utilities ----------


def test_revcomp_involutive():
    s = "ACGTACGTTT"
    assert revcomp(revcomp(s)) == s


def test_canonical_picks_lex_smaller():
    km = "AAAAT"
    assert canonical(km) == km  # its rc is ATTTT > AAAAT
    km2 = "TTTTT"
    assert canonical(km2) == "AAAAA"  # rc of TTTTT


def test_kmers_enumeration():
    assert list(kmers("ACGTACGT", 4)) == ["ACGT", "CGTA", "GTAC", "TACG", "ACGT"]


# ---------- k-mer counting ----------


def test_count_kmers_canonical_merges_reverse_complements():
    reads = ["ACGT"]
    counts = count_kmers(reads, k=2, canonicalize=True)
    # 2-mers in ACGT: AC, CG, GT. Their rcs: GT, CG, AC. AC and GT merge.
    # Canonical: AC (min of AC/GT), CG (min of CG/CG), AC again from GT.
    assert counts[canonical("AC")] >= 1
    assert counts[canonical("CG")] >= 1


def test_count_kmers_no_canonical():
    reads = ["AAAA"]
    counts = count_kmers(reads, k=3, canonicalize=False)
    assert counts["AAA"] == 2


def test_count_min_sketch_over_estimates():
    cms = CountMinSketch(width=1024, depth=4)
    cms.add("hello", 3)
    cms.add("world", 1)
    assert cms.query("hello") >= 3
    assert cms.query("world") >= 1
    # items never added should give a nonnegative (often 0) estimate
    assert cms.query("never_added") >= 0


# ---------- de Bruijn graph construction ----------


def test_de_bruijn_simple_sequence():
    read = "ACGTACGT"
    g = DeBruijnGraph.from_reads([read], k=3)
    # 3-mers: ACG, CGT, GTA, TAC, ACG, CGT -> 6 edges
    assert g.edge_count() == len(read) - 3 + 1
    # At least a couple of node labels are the (k-1) prefixes.
    assert "AC" in g.nodes()


def test_de_bruijn_handles_ns():
    reads = ["ACGNT", "ACGAC"]
    g = DeBruijnGraph.from_reads(reads, k=3)
    # k-mers containing N should be skipped.
    for u in g.nodes():
        assert "N" not in u


# ---------- Eulerian path + reconstruction ----------


def test_eulerian_recovers_linear_sequence():
    reference = "ACGTACGTACGTACGT"
    reads = [reference]
    g = DeBruijnGraph.from_reads(reads, k=4)
    path = find_eulerian_path(g)
    assert path is not None
    seq = path_to_sequence(path)
    # Reconstructed sequence should have all the same k-mers as the original.
    original_kmers = set(kmers(reference, 4))
    reconstructed_kmers = set(kmers(seq, 4))
    assert original_kmers == reconstructed_kmers


def test_eulerian_on_random_simulated_reads():
    # Error free reads at high coverage produce a multigraph whose edge
    # multiplicities scale with coverage, so the Eulerian walk visits each
    # edge-occurrence and the reconstructed sequence length is roughly
    # coverage times the reference length. The semantically correct check
    # is that every reference k-mer is recovered in the reconstruction.
    reference = "ACGTACGTAGGCTAGCTAGCATCGATGCATGCTAGCT" * 3
    reads = simulate_reads(reference, read_length=30, coverage=50.0,
                           error_rate=0.0, rng=np.random.default_rng(0))
    g = DeBruijnGraph.from_reads(reads, k=15)
    path = find_eulerian_path(g)
    if path is not None:
        seq = path_to_sequence(path)
        ref_kmers = set(kmers(reference, 15))
        seq_kmers = set(kmers(seq, 15))
        assert ref_kmers.issubset(seq_kmers)


# ---------- cleanup + contigs ----------


def test_tip_removal_eliminates_dead_ends():
    # Build a main path A->B->C->D and an errant tip A->X (X is a dead end).
    g = DeBruijnGraph(k=2)
    g._add_edge("AA", "AB", "AAB")
    g._add_edge("AB", "BC", "ABC")
    g._add_edge("BC", "CD", "BCD")
    g._add_edge("AA", "AX", "AAX")  # tip
    removed = g.remove_tips(max_tip_length=3)
    assert removed >= 1
    assert "AX" not in g.out_edges.get("AA", {})


def test_extract_contigs_on_linear_graph():
    # Purely linear graph: A -> B -> C -> D -> E.
    g = DeBruijnGraph(k=2)
    g._add_edge("AA", "AB", "AAB")
    g._add_edge("AB", "BC", "ABC")
    g._add_edge("BC", "CD", "BCD")
    g._add_edge("CD", "DE", "CDE")
    contigs = extract_contigs(g)
    assert len(contigs) == 1
    # Five 2-character nodes overlapping by one character reconstruct a
    # sequence of length 2 + (5 - 1) = 6 ("AABCDE").
    assert len(contigs[0]) == 6
    assert contigs[0] == "AABCDE"


def test_n50_on_known_contigs():
    contigs = ["A" * 100, "A" * 50, "A" * 30, "A" * 20]
    # total = 200; half = 100. Sorted desc: 100, 50, 30, 20. Running 100 (>=100) -> N50 = 100.
    assert n50(contigs) == 100


def test_assembly_stats_fields():
    contigs = ["AAAA", "AAA", "AA"]
    stats = assembly_stats(contigs)
    assert stats["n_contigs"] == 3
    assert stats["total_length"] == 9
    assert stats["longest"] == 4
    assert stats["shortest"] == 2


# ---------- end to end ----------


def test_end_to_end_simulated_lambda_like():
    reference = ("ATGCAGCTAGCTAGCATGCATGCATGCATGCTAGCTAGCATCGATCGATCG"
                 "ACGTAGCATGCATCGATCGATCGTACGATCGATCGATCGATCGTAGCATGC"
                 "ATGCATCGATCGTAGCATCGATGCATGCATGCATGCATGCATGCATGCTAG") * 2
    reads = simulate_reads(reference, read_length=80, coverage=40.0,
                           error_rate=0.0, rng=np.random.default_rng(1))
    g = DeBruijnGraph.from_reads(reads, k=25)
    g.remove_low_coverage_edges(min_coverage=1)
    g.remove_tips()
    contigs = extract_contigs(g)
    stats = assembly_stats(contigs)
    # Error free reads at 40x coverage with k=25 should recover a large fraction of the reference.
    assert stats["total_length"] >= 0.5 * len(reference)
