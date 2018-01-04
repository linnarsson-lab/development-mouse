from typing import *
import os
import logging
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
import networkx as nx
import cytograph as cg
import loompy
from matplotlib.colors import LinearSegmentedColormap
import numpy_groupies.aggregate_numpy as npg
import scipy.cluster.hierarchy as hc
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as path_effects
import matplotlib.colors as mcolors
from matplotlib.colors import colorConverter
from matplotlib.collections import LineCollection
from sklearn.neighbors import BallTree, NearestNeighbors, kneighbors_graph


def plot_abs_graph(ds: loompy.LoomConnection, dsagg: loompy.LoomConnection, out_file: str=None, tags: List[str] = None) -> None:
    logging.info("Loading graph")
    n_cells = ds.shape[1]
    cells = np.where(ds.col_attrs["_Valid"] == 1)[0]
    pos = np.vstack((ds.col_attrs["_X"], ds.col_attrs["_Y"])).transpose()
    labels = ds.col_attrs["Clusters"]

    # Compute a good size for the markers, based on local density
    logging.info("Computing node size")
    min_pts = 50
    eps_pct = 60
    nn = NearestNeighbors(n_neighbors=min_pts, algorithm="ball_tree", n_jobs=4)
    nn.fit(pos)
    knn = nn.kneighbors_graph(mode='distance')
    k_radius = knn.max(axis=1).toarray()
    epsilon = 24 * np.percentile(k_radius, eps_pct)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    # Draw (faded) single cell nodes
    logging.info("Drawing nodes")
    colors20 = np.vstack((plt.cm.Vega20b(np.linspace(0., 1, 20))[::2], plt.cm.Vega20c(np.linspace(0, 1, 20))[1::2]))
    plots = []
    names = []
    for i in range(max(labels) + 1):
        cluster = labels == i
        plots.append(plt.scatter(x=pos[cluster, 0], y=pos[cluster, 1], c=colors20[np.mod(i, 20)], marker='.', lw=0, s=epsilon, alpha=0.05, zorder=-1))
        if tags is not None:
            names.append(str(i) + " " + tags[i].replace("\n", " "))
        else:
            names.append(str(i))

    # Retrieve graph info from AggregateL1 results
    (a, b, w) = dsagg.get_edges("AbstractGraph", axis=1)

    # Coordinates of the abstracted nodes
    pos_ag = np.zeros((max(labels) + 1, 2))
    for lbl in range(0, max(labels) + 1):
        pos_ag[lbl, :] = np.median(pos[np.where(labels == lbl)[0]], axis=0)

    # Draw abstracted graph edges
    logging.info("Drawing edges")
    confidence_bins_pars = {(0, 0.25): ('dotted', 0.5, 0.3),
                            (0.25, 0.5): ('dotted', 1, 0.3),
                            (0.5, 0.75): ('dashed', 2, 0.5),
                            (0.75, 0.9): ('solid', 3, 0.7),
                            (0.9, 1): ('solid', 7, 0.9)}
    for (m, M), (line_style, line_width, alpha) in confidence_bins_pars.items():
        is_in_bin = (w > m) & (w <= M)
        if np.any(is_in_bin):
            lc = LineCollection(zip(pos_ag[a[is_in_bin]], pos_ag[b[is_in_bin]]), linestyles=line_style, linewidths=line_width, zorder=10000, color='grey', alpha=alpha)
            ax.add_collection(lc)

    # Draw abstracted graph nodes and their IDS
    logging.info("Drawing cluster graph nodes and IDs")
    plt.scatter(pos_ag[:, 0], pos_ag[:, 1], c=colors20[np.mod(np.arange(pos_ag.shape[0]), 20)], marker='o', s=150, alpha=0.8, zorder=20000)
    for lbl in range(0, max(labels) + 1):
        (x, y) = pos_ag[lbl, :]  # np.median(pos[np.where(labels == lbl)[0]], axis=0)
        ax.text(x, y, str(lbl), fontsize=12, bbox=dict(facecolor='white', alpha=0.5, ec='none'))

    # Drawing legend
    logging.info("Drawing legend")
    plt.legend(plots, names, scatterpoints=1, markerscale=2, loc='upper left', bbox_to_anchor=(1, 1), fancybox=True, framealpha=0.5, fontsize=10)
    
    # Save to file
    if out_file is not None:
        logging.info("Saving to file")
        fig.savefig(out_file, format="png", dpi=144, bbox_inches='tight')
        plt.close()


def plot_velocity_summary(vlm: Any, confidence: np.ndarray, significant: np.ndarray, trans: np.ndarray,
                          expected_tr: np.ndarray, out_file: str=None, tags: List[str] = None) -> None:
    logging.info("Loading info from VelocytoLoom object")
    n_cells = vlm.S.shape[1]
    cells = np.ones(n_cells, dtype=bool)
    labels = vlm.cluster_ix
    pos = vlm.ts

    # Compute a good size for the markers, based on local density
    logging.info("Computing node size")
    min_pts = 50
    eps_pct = 60
    nn = NearestNeighbors(n_neighbors=min_pts, algorithm="ball_tree", n_jobs=4)
    nn.fit(pos)
    knn = nn.kneighbors_graph(mode='distance')
    k_radius = knn.max(axis=1).toarray()
    epsilon = 24 * np.percentile(k_radius, eps_pct)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    # Draw (faded) single cell nodes
    logging.info("Drawing nodes")
    colors20 = np.vstack((plt.cm.Vega20b(np.linspace(0., 1, 20))[::2], plt.cm.Vega20c(np.linspace(0, 1, 20))[1::2]))
    plots = []
    names = []
    for i in range(max(labels) + 1):
        cluster = labels == i
        plots.append(plt.scatter(x=pos[cluster, 0], y=pos[cluster, 1], c=colors20[np.mod(i, 20)], marker='.', lw=0, s=epsilon, alpha=0.05, zorder=-1))
        if tags is not None:
            names.append(str(i) + " " + tags[i].replace("\n", " "))
        else:
            names.append(str(i))

    # Retrieve graph info from AggregateL1 results
    # (a, b, w) = dsagg.get_edges("AbstractGraph", axis=1)
    confidence_ = np.copy(confidence)
    confidence_[confidence_ < 0.2] = 0  # NOTE: this threshold is a bit arbitrary
    abs_graph = sparse.coo_matrix(confidence_)
    a, b, w = abs_graph.row, abs_graph.col, abs_graph.data

    # Coordinates of the abstracted nodes
    pos_ag = np.zeros((max(labels) + 1, 2))
    for lbl in range(0, max(labels) + 1):
        pos_ag[lbl, :] = np.median(pos[np.where(labels == lbl)[0]], axis=0)

    # Draw abstracted graph edges
    logging.info("Drawing edges")
    confidence_bins_pars = {(0, 0.25): ('dotted', 0.5, 0.3),
                            (0.25, 0.5): ('dotted', 1, 0.3),
                            (0.5, 0.75): ('dashed', 2, 0.5),
                            (0.75, 0.9): ('solid', 3, 0.7),
                            (0.9, 1): ('solid', 7, 0.9)}
    for (m, M), (line_style, line_width, alpha) in confidence_bins_pars.items():
        is_in_bin = (w > m) & (w <= M)
        if np.any(is_in_bin):
            lc = LineCollection(zip(pos_ag[a[is_in_bin]], pos_ag[b[is_in_bin]]), linestyles=line_style, linewidths=line_width, zorder=10000, color='grey', alpha=alpha)
            ax.add_collection(lc)

    # Draw arrows for transition
    c, d = np.where(significant)
    for i in range(len(c)):
        lbl_from, lbl_to = c[i], d[i]
        (x, y) = pos_ag[lbl_from, :]
        (x2, y2) = pos_ag[lbl_to, :]
        dx = x2 - x
        dy = y2 - y
        plt.arrow(x, y, dx, dy, zorder=200000, length_includes_head=True, width=1)

    # Draw abstracted graph nodes and their IDS
    logging.info("Drawing cluster graph nodes and IDs")
    plt.scatter(pos_ag[:, 0], pos_ag[:, 1], c=colors20[np.mod(np.arange(pos_ag.shape[0]), 20)], marker='o', s=250, alpha=0.8, zorder=20000)
    for lbl in range(0, max(labels) + 1):
        (x, y) = pos_ag[lbl, :]  # np.median(pos[np.where(labels == lbl)[0]], axis=0)
        ax.text(x, y, str(lbl), fontsize=12, bbox=dict(facecolor='white', alpha=0.5, ec='none'))

    # Drawing lagend
    logging.info("Drawing legend")
    plt.legend(plots, names, scatterpoints=1, markerscale=2, loc='upper left', bbox_to_anchor=(1, 1), fancybox=True, framealpha=0.5, fontsize=10)
    
    # Save to file
    if out_file is not None:
        logging.info("Saving to file")
        fig.savefig(out_file, format="png", dpi=144, bbox_inches='tight')
        plt.close()
