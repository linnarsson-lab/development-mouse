from typing import *
import os
import scanpy.api as sc
import luigi
import development_mouse as dm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cytograph as cg
import logging
import loompy
import numpy as np

class PlotAGA(luigi.Task):
    """
    Luigi task to plot scanpy AGA
    """
    tissue = luigi.Parameter(description="name of the tissue from tool specification file.",
                             always_in_help=True)

    def requires(self) -> List[luigi.Task]:
        return [dm.ClusterL1(tissue=self.tissue),
		dm.RunAGA(tissue=self.tissue)]

    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, f"AGA_{self.tissue}.png"))

    def run(self) -> None:
        logging.info("Plotting AGA")
        with self.output().temporary_path() as out_file:
                adata = sc.read(self.input()[1].fn)
                ds = loompy.connect(self.input()[0].fn)
                label = ds.col_attrs["Clusters"]
                clusters = np.unique(label)
                cMap = colors.ListedColormap(plt.cm.tab20(clusters))

                fig = plt.figure(figsize=(10, 5))
                ax1 = fig.add_subplot(121)
                ax2 = fig.add_subplot(122)
                ax1.scatter(ds.col_attrs["_X"], ds.col_attrs["_Y"], c=label, cmap=cMap, s=2)
                for i in clusters:
                    x = np.median(ds.col_attrs["_X"][label == i])
                    y = np.median(ds.col_attrs["_Y"][label == i])
                    ax1.text(x, y, np.unique(label)[i], fontsize=10, bbox={"facecolor":"w", "alpha":0.6})
                ax1.set_title('cytograph clusters')
                sc.pl.aga_graph(adata, ax=ax2, node_size_scale=2, edge_width_scale=1.5, fontsize=15)
                fig.tight_layout()
                fig.savefig(out_file, format='png', dpi=300)           
