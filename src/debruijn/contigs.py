"""Contig extraction from a (cleaned) de Bruijn graph.

A contig is a maximal non branching path: a walk along the graph that
starts at a branching node (or dead end), follows a chain of nodes each
with in degree 1 and out degree 1, and ends at the next branching node
or dead end. Every edge in the graph belongs to exactly one such path.

This is the classical "unitig graph" construction and is how assemblers
like Velvet output their contigs.
"""

from __future__ import annotations

from typing import List, Set, Tuple

from .eulerian import path_to_sequence
from .graph import DeBruijnGraph


def _is_branching_or_endpoint(g: DeBruijnGraph, u: str) -> bool:
    return g.out_degree(u) != 1 or g.in_degree(u) != 1


def extract_contigs(g: DeBruijnGraph) -> List[str]:
    """Return a list of contig sequences obtained by collapsing all
    maximal non branching paths in the graph."""
    visited_edges: Set[Tuple[str, str]] = set()
    contigs: List[str] = []

    # Walk out from every branching or endpoint node along each outgoing edge.
    for u in list(g.nodes()):
        if not _is_branching_or_endpoint(g, u):
            continue
        for v in list(g.out_edges.get(u, {}).keys()):
            if (u, v) in visited_edges:
                continue
            path = [u, v]
            visited_edges.add((u, v))
            cur = v
            while (
                g.out_degree(cur) == 1
                and g.in_degree(cur) == 1
                and not _is_branching_or_endpoint(g, cur)
            ):
                nxts = list(g.out_edges[cur].keys())
                if not nxts:
                    break
                nxt = nxts[0]
                if (cur, nxt) in visited_edges:
                    break
                visited_edges.add((cur, nxt))
                path.append(nxt)
                cur = nxt
            contigs.append(path_to_sequence(path))

    # Isolated cycles with no branching nodes — walk from any remaining node.
    for u in list(g.nodes()):
        for v in list(g.out_edges.get(u, {}).keys()):
            if (u, v) in visited_edges:
                continue
            path = [u, v]
            visited_edges.add((u, v))
            cur = v
            while cur != u:
                nxts = list(g.out_edges.get(cur, {}).keys())
                if not nxts:
                    break
                nxt = nxts[0]
                if (cur, nxt) in visited_edges:
                    break
                visited_edges.add((cur, nxt))
                path.append(nxt)
                cur = nxt
            contigs.append(path_to_sequence(path))

    # Dedup and sort by length
    contigs = sorted(set(contigs), key=len, reverse=True)
    return contigs


def n50(contigs: List[str]) -> int:
    """N50: the length N such that contigs of length >= N sum to at
    least half the total assembly length."""
    if not contigs:
        return 0
    sorted_lens = sorted((len(c) for c in contigs), reverse=True)
    total = sum(sorted_lens)
    half = total / 2
    running = 0
    for L in sorted_lens:
        running += L
        if running >= half:
            return L
    return sorted_lens[-1]


def assembly_stats(contigs: List[str]) -> dict:
    lens = sorted((len(c) for c in contigs), reverse=True)
    return {
        "n_contigs": len(contigs),
        "total_length": sum(lens),
        "longest": lens[0] if lens else 0,
        "shortest": lens[-1] if lens else 0,
        "n50": n50(contigs),
        "mean_length": (sum(lens) / len(lens)) if lens else 0.0,
    }
