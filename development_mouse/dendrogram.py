import numpy as np
import loompy
import scipy.cluster.hierarchy as hc
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from polo import optimal_leaf_ordering


def plot_dendrogram(dsagg: loompy.LoomConnection) -> None:
    selection = dsagg.row_attrs["_Selected"] == 1
    data = np.log(dsagg[:, :] + 1)[selection, :].T
    D = pdist(data, 'euclidean')
    Z = hc.linkage(D, 'ward')
    opt_Z = optimal_leaf_ordering(Z, D)
    opt_order = hc.leaves_list(opt_Z)

    # Plotting
    ax = plt.subplot(111)
    hc.dendrogram(opt_Z, ax=ax, link_color_func=lambda k: 'k')
    ax.set_xticklabels(data[opt_order].reshape(-1))
    ax.set_xticks([])
    ax.set_yticks([])
