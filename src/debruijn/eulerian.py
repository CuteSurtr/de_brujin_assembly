"""Hierholzer's algorithm for finding Eulerian circuits and paths in a
de Bruijn graph. Runs in O(E) time.

An Eulerian path visits every edge exactly once. A directed multigraph
has an Eulerian circuit iff every node has equal in and out degree and
the graph (restricted to nodes with at least one edge) is connected. It
has an Eulerian path iff exactly one node has out_degree - in_degree = 1
(the start) and exactly one has in_degree - out_degree = 1 (the end),
with all other nodes balanced.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from .graph import DeBruijnGraph


def find_eulerian_path(graph: DeBruijnGraph) -> Optional[List[str]]:
    """Return a node sequence representing an Eulerian path through the
    graph, or None if none exists. Multi edges between the same pair of
    nodes are consumed one at a time; since our DeBruijnGraph stores
    edge multiplicities as weights, we decrement the weight each time an
    edge is traversed and remove the edge when weight reaches zero.

    In a directed multigraph, the degree counts used to determine the
    start must include multiplicities. A node has ``total_out_degree``
    equal to the sum of outgoing edge weights and similarly for the
    incoming side. An Eulerian path exists iff exactly one node has
    total_out_degree - total_in_degree = +1 (the start) and exactly one
    has -1 (the end), with all other nodes balanced; or all nodes are
    balanced in which case there is an Eulerian circuit.
    """

    def total_out(u: str) -> int:
        return sum(w for _, w in graph.out_edges.get(u, {}).values())

    def total_in(u: str) -> int:
        return sum(w for _, w in graph.in_edges.get(u, {}).values())

    # Prefer a node with an excess of one outgoing edge; this is the
    # unique Eulerian-path start when it exists.
    start = None
    for u in graph.nodes():
        if total_out(u) - total_in(u) == 1:
            start = u
            break
    if start is None:
        for u in graph.nodes():
            if total_out(u) > 0:
                start = u
                break
    if start is None:
        return None

    # Make a mutable copy of edge weights.
    out_weights: Dict[str, Dict[str, int]] = {
        u: {v: w for v, (_, w) in nbrs.items()}
        for u, nbrs in graph.out_edges.items()
    }

    path: List[str] = []
    stack: List[str] = [start]
    while stack:
        u = stack[-1]
        if u in out_weights and any(w > 0 for w in out_weights[u].values()):
            # Pick any edge with remaining weight.
            for v, w in out_weights[u].items():
                if w > 0:
                    out_weights[u][v] -= 1
                    stack.append(v)
                    break
        else:
            path.append(stack.pop())

    path.reverse()
    # Check all edges consumed.
    total_remaining = sum(sum(d.values()) for d in out_weights.values())
    if total_remaining > 0:
        return None
    return path


def path_to_sequence(path: List[str]) -> str:
    """Convert a path of overlapping (k-1)-mers into the concatenated
    sequence they represent. Each successive node adds exactly one
    character."""
    if not path:
        return ""
    out = [path[0]]
    for i in range(1, len(path)):
        out.append(path[i][-1])
    return "".join(out)
