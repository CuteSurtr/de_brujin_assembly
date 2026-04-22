"""De Bruijn graph construction and cleanup for genome assembly.

Convention: a de Bruijn graph of order k has
  * nodes = distinct (k-1)-mers observed in the reads,
  * edges = distinct k-mers observed in the reads, connecting their
    (k-1)-prefix to their (k-1)-suffix,
  * edge weights = k-mer coverage (number of reads contributing).

We store the graph as two dicts: ``out_edges[node] -> {neighbor: (edge_label, weight)}``
and ``in_edges[node] -> {predecessor: (edge_label, weight)}``. The edge
label is the full k-mer string so that contig extraction can reconstruct
sequences exactly.

For double stranded DNA, we canonicalize k-mers (lexicographic min of
k-mer and reverse complement), which corresponds to treating both
strands as contributing to the same edge. Traversal then emits the
canonical path; reconstructing the forward strand is a separate post
processing step.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Set, Tuple

from .utils import canonical, kmers


@dataclass
class DeBruijnGraph:
    k: int
    out_edges: Dict[str, Dict[str, Tuple[str, int]]] = field(default_factory=lambda: defaultdict(dict))
    in_edges: Dict[str, Dict[str, Tuple[str, int]]] = field(default_factory=lambda: defaultdict(dict))

    @classmethod
    def from_reads(cls, reads: Iterable[str], k: int, canonicalize: bool = False) -> "DeBruijnGraph":
        """Build a de Bruijn graph of order k from a read set. By default
        we do NOT canonicalize (so the graph respects strand orientation
        as seen in the reads); pass ``canonicalize=True`` for a canonical
        graph that merges both strands."""
        g = cls(k=k)
        for read in reads:
            for km in kmers(read, k):
                if canonicalize:
                    km = canonical(km)
                u, v = km[:-1], km[1:]
                g._add_edge(u, v, km)
        return g

    def _add_edge(self, u: str, v: str, label: str) -> None:
        if v in self.out_edges[u]:
            old_label, w = self.out_edges[u][v]
            self.out_edges[u][v] = (old_label, w + 1)
            self.in_edges[v][u] = (old_label, w + 1)
        else:
            self.out_edges[u][v] = (label, 1)
            self.in_edges[v][u] = (label, 1)
        # Ensure the reverse direction dict exists even if no outgoing edges
        self.out_edges.setdefault(v, {})
        self.in_edges.setdefault(u, {})

    def _remove_edge(self, u: str, v: str) -> None:
        if v in self.out_edges.get(u, {}):
            del self.out_edges[u][v]
        if u in self.in_edges.get(v, {}):
            del self.in_edges[v][u]

    def nodes(self) -> Set[str]:
        return set(self.out_edges.keys()) | set(self.in_edges.keys())

    def out_degree(self, u: str) -> int:
        return len(self.out_edges.get(u, {}))

    def in_degree(self, u: str) -> int:
        return len(self.in_edges.get(u, {}))

    def edge_count(self) -> int:
        """Total edge count including multiplicity (i.e., sum of edge weights)."""
        return sum(w for nbrs in self.out_edges.values() for (_, w) in nbrs.values())

    def distinct_edge_count(self) -> int:
        """Number of distinct (u, v) pairs, ignoring multiplicity."""
        return sum(len(d) for d in self.out_edges.values())

    # ---------- cleanup ----------

    def remove_low_coverage_edges(self, min_coverage: int) -> int:
        """Drop any edge whose weight is below ``min_coverage``. Returns
        the number of edges removed."""
        to_remove: List[Tuple[str, str]] = []
        for u, nbrs in list(self.out_edges.items()):
            for v, (_, w) in list(nbrs.items()):
                if w < min_coverage:
                    to_remove.append((u, v))
        for u, v in to_remove:
            self._remove_edge(u, v)
        return len(to_remove)

    def remove_tips(self, max_tip_length: Optional[int] = None) -> int:
        """Iteratively prune tips shorter than ``max_tip_length`` edges.

        A forward tip is a chain ending in a node with out_degree 0 whose
        walk back terminates at a branch point (a node with in-degree or
        out-degree greater than one). A reverse tip is the symmetric case
        starting from a node with in_degree 0. If the tip has strictly
        fewer edges than ``max_tip_length`` (Velvet 2008 uses 2k, with
        strict inequality), it is removed.

        Returns the number of edges removed.
        """
        if max_tip_length is None:
            max_tip_length = 2 * self.k
        removed = 0
        changed = True
        while changed:
            changed = False
            # Forward tips: walk backward from out-degree-zero nodes.
            for dead in [n for n in list(self.nodes())
                         if self.out_degree(n) == 0 and self.in_degree(n) >= 1]:
                chain = [dead]
                cur = dead
                while len(chain) <= max_tip_length:
                    preds = list(self.in_edges.get(cur, {}).keys())
                    if len(preds) != 1:
                        break
                    p = preds[0]
                    if self.out_degree(p) > 1 or self.in_degree(p) > 1:
                        chain.append(p)
                        break
                    chain.append(p)
                    cur = p
                last = chain[-1]
                at_branch = self.out_degree(last) > 1 or self.in_degree(last) > 1
                tip_edges = len(chain) - 1
                if at_branch and 1 <= tip_edges < max_tip_length:
                    for i in range(len(chain) - 1):
                        self._remove_edge(chain[i + 1], chain[i])
                    removed += tip_edges
                    changed = True
            # Reverse tips: walk forward from in-degree-zero nodes.
            for dead in [n for n in list(self.nodes())
                         if self.in_degree(n) == 0 and self.out_degree(n) >= 1]:
                chain = [dead]
                cur = dead
                while len(chain) <= max_tip_length:
                    succs = list(self.out_edges.get(cur, {}).keys())
                    if len(succs) != 1:
                        break
                    s = succs[0]
                    if self.out_degree(s) > 1 or self.in_degree(s) > 1:
                        chain.append(s)
                        break
                    chain.append(s)
                    cur = s
                last = chain[-1]
                at_branch = self.out_degree(last) > 1 or self.in_degree(last) > 1
                tip_edges = len(chain) - 1
                if at_branch and 1 <= tip_edges < max_tip_length:
                    for i in range(len(chain) - 1):
                        self._remove_edge(chain[i], chain[i + 1])
                    removed += tip_edges
                    changed = True
        return removed

    def pop_bubbles(self, max_coverage_ratio: float = 0.5) -> int:
        """Identify pairs of parallel paths between two shared nodes and
        remove the lower coverage path if its coverage ratio to the
        higher is below ``max_coverage_ratio``."""
        removed = 0
        changed = True
        while changed:
            changed = False
            for u in list(self.out_edges.keys()):
                neighbors = list(self.out_edges[u].items())
                if len(neighbors) < 2:
                    continue
                # Look for simple bubbles: two short paths from u that
                # converge at the same node w.
                for i in range(len(neighbors)):
                    for j in range(i + 1, len(neighbors)):
                        (v1, (_, w1)) = neighbors[i]
                        (v2, (_, w2)) = neighbors[j]
                        nxts1 = list(self.out_edges.get(v1, {}).keys())
                        nxts2 = list(self.out_edges.get(v2, {}).keys())
                        if (
                            len(nxts1) == 1 and len(nxts2) == 1 and nxts1[0] == nxts2[0]
                            and self.in_degree(v1) == 1 and self.in_degree(v2) == 1
                            and self.out_degree(v1) == 1 and self.out_degree(v2) == 1
                        ):
                            if w1 >= w2 and w2 / max(w1, 1) <= max_coverage_ratio:
                                self._remove_edge(u, v2)
                                self._remove_edge(v2, nxts2[0])
                                removed += 2
                                changed = True
                                break
                            elif w2 > w1 and w1 / max(w2, 1) <= max_coverage_ratio:
                                self._remove_edge(u, v1)
                                self._remove_edge(v1, nxts1[0])
                                removed += 2
                                changed = True
                                break
                    if changed:
                        break
                if changed:
                    break
        return removed
