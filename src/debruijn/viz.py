from __future__ import annotations

from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np


def plot_kmer_histogram(h: np.ndarray, ax: Optional[plt.Axes] = None, title: str = "k-mer frequency histogram"):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.bar(np.arange(len(h)), h, color="tab:blue", width=1.0)
    ax.set_xlabel("observed k-mer count")
    ax.set_ylabel("distinct k-mers with this count")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    return ax


def plot_contig_lengths(contig_lens: Sequence[int], ax: Optional[plt.Axes] = None, title: str = "Contig length distribution"):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 3.5))
    lens_sorted = sorted(contig_lens, reverse=True)
    ax.plot(range(1, len(lens_sorted) + 1), lens_sorted, "o-", markersize=4, color="tab:green")
    ax.set_xlabel("contig rank")
    ax.set_ylabel("contig length (bp)")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    return ax


def plot_graph_degree_distribution(graph, ax: Optional[plt.Axes] = None, title: str = "de Bruijn graph degree distribution"):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 3.5))
    in_degs = [graph.in_degree(u) for u in graph.nodes()]
    out_degs = [graph.out_degree(u) for u in graph.nodes()]
    max_d = max(max(in_degs, default=0), max(out_degs, default=0)) + 1
    bins = np.arange(max_d + 1)
    ax.hist([in_degs, out_degs], bins=bins, label=["in-degree", "out-degree"], color=["tab:blue", "tab:orange"])
    ax.set_xlabel("degree")
    ax.set_ylabel("number of nodes")
    ax.set_title(title)
    ax.legend()
    return ax
