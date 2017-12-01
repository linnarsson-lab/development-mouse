from typing import *
import os
import loompy
from scipy import sparse
import numpy as np
import networkx as nx
import cytograph as cg
import development_mouse as dm
import luigi
import matplotlib.colors as colors
import matplotlib.pyplot as plt


class ExportL1(luigi.Task):
    """Luigi Task to autoannotate and export summary files
    """
    tissue = luigi.Parameter(description="name of the tissue from tool specification file.",
                             always_in_help=True)

    n_markers = luigi.IntParameter(default=10,
                                   description="number of markers visualized in the heatmap")

    def requires(self) -> Dict[str, luigi.Task]:
        """
        Arguments
        ---------
        `AggregateL1`:
            passing ``tissue``
        `ClusterL1`:
            passing ``tissue``
        """
        return {f"AggregateL1(tissue={self.tissue})": dm.AggregateL1(tissue=self.tissue),
                f"ClusterL1(tissue={self.tissue})": dm.ClusterL1(tissue=self.tissue),
                "NameQualityClusters": dm.NameQualityClusters()}

    def output(self) -> luigi.Target:
        """
        Returns
        -------
        folder: ``L1_[TISSUE]_exported``:
            Note this is kind of a hack to luigi, single files will not be regenerated but whole folder will.
        """
        return luigi.LocalTarget(os.path.join(dm.paths().build, "L1_" + self.tissue + "_exported"))

    def run(self) -> None:
        """
        Reads the output of `AggregateL1` and does:
            - runs the `cytograph.AutoAnnotator`
            - exports  ``L1_[TISSUE]_expression.tab``, ``L1_[TISSUE]_enrichment.tab``, ``L1_[TISSUE]_trinaries.tab``
            - uses `cytograph.plot_graph` to plot ``L1_[TISSUE]_manifold.aa.png``, ``L1_[TISSUE]_manifold.aaa.png``
            - uses `cytograph.plot_markerheatmap` to plot ``L1_[TISSUE]_heatmap.pdf``

        Other Parameters
        ----------------
        `memory_config.memory.batch`: int
        the number of columns/rows used by batchscan

        """
        logging = cg.logging(self, True)
        logging.info("Exporting cluster data")
        with self.output().temporary_path() as out_dir:
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            dsagg = loompy.connect(self.input()[f"AggregateL1(tissue={self.tissue})"].fn)

            dsagg.export(os.path.join(out_dir, "L1_" + self.tissue + "_expression.tab"))
            dsagg.export(os.path.join(out_dir, "L1_" + self.tissue + "_enrichment.tab"), layer="enrichment")
            dsagg.export(os.path.join(out_dir, "L1_" + self.tissue + "_enrichment_q.tab"), layer="enrichment_q")
            dsagg.export(os.path.join(out_dir, "L1_" + self.tissue + "_trinaries.tab"), layer="trinaries")

            ds = loompy.connect(self.input()[f"ClusterL1(tissue={self.tissue})"].fn)

            logging.info("Plotting manifold graph with auto-annotation")
            tags = list(dsagg.col_attrs["AutoAnnotation"])
            cg.plot_graph(ds, os.path.join(out_dir, "L1_" + self.tissue + "_manifold.aa.png"), tags)
            cg.plot_graph_age(ds, os.path.join(out_dir, "L1_" + self.tissue + "_manifold.age.png"), tags)

            logging.info("Plotting manifold graph with auto-auto-annotation")
            tags = list(dsagg.col_attrs["MarkerGenes"])
            cg.plot_graph(ds, os.path.join(out_dir, "L1_" + self.tissue + "_manifold.aaa.png"), tags)

            logging.info("Plotting marker heatmap")
            cg.plot_markerheatmap(ds, dsagg, n_markers_per_cluster=self.n_markers, out_file=os.path.join(out_dir, "L1_" + self.tissue + "_heatmap.pdf"))

            logging.info("Plotting quality class on t-SNE")
            tags = list(dsagg.col_attrs["AutoAnnotation"])
            cluster_mapping = {int(i.split(":")[0]): i.split(":")[1] 
                                    for i in open(self.input()["NameQualityClusters"].fn).read().rstrip().split()}
            n_cells = ds.shape[1]
            valid = ds.col_attrs["_Valid"].astype('bool')
            
            (a, b, w) = ds.get_edges("MKNN", axis=1)
            mknn = sparse.coo_matrix((w, (a, b)), shape=(n_cells, n_cells)).tocsr()[valid, :][:, valid]
            sfdp = np.vstack((ds.col_attrs["_X"], ds.col_attrs["_Y"])).transpose()[valid, :]
            orderx = np.argsort(sfdp[:, 0], kind="mergesort")
            ordery = np.argsort(sfdp[:, 1], kind="mergesort")
            orderfin = orderx[ordery]
            sfdp_original = sfdp.copy()  
            sfdp = sfdp[orderfin, :]
            labels = ds.col_attrs["Clusters"][valid][orderfin]
            quality = ds.col_attrs["QualityClass"].astype(float)[valid][orderfin]
            
            fig = plt.figure(figsize=(10, 10))
            g = nx.from_scipy_sparse_matrix(mknn)
            ax = fig.add_subplot(111)

            nx.draw_networkx_edges(g, pos=sfdp_original, alpha=0.1, width=0.1, edge_color='gray')
            block_colors = plt.cm.tab20(quality / np.max(quality))
            nx.draw_networkx_nodes(g, pos=sfdp, node_color=block_colors, node_size=10, alpha=0.6, linewidths=0)

            mg_pos = []
            for lbl in range(0, max(labels) + 1):
                if np.sum(labels == lbl) == 0:
                    continue
                (x, y) = np.median(sfdp[np.where(labels == lbl)[0]], axis=0)
                mg_pos.append((x, y))
                text = "#" + str(lbl)
                if len(tags[lbl]) > 0:
                    text += "\n" + tags[lbl]
                ax.text(x, y, text, fontsize=8, weight="bold", bbox=dict(facecolor='gray', alpha=0.3, ec='none'))
            ax.axis('off')
            levels = np.unique(quality)
            for il, lev in enumerate(levels):
                ax.add_patch(
                    plt.Rectangle(
                        (0.90, 0.7 + il * 0.016), 0.014, 0.014,
                        color=plt.cm.tab20(lev / np.max(levels)),
                        clip_on=0, transform=ax.transAxes))
                ax.text(0.93, 0.703 + il * 0.016, cluster_mapping[lev], transform=ax.transAxes)
            plt.tight_layout()
            fig.savefig(os.path.join(out_dir, "L1_" + self.tissue + "_quality_tsne.png"), format="png", dpi=300)
            plt.close()
                    
            logging.info("Plotting quality class in pie chart")
            plt.figure(None, (10, 10))
            labels = ds.col_attrs["QualityClass"].astype(int)
            unique, counts = np.unique(labels, return_counts=True)
            labelnames = [cluster_mapping[ix] for ix in unique]
            patches, texts = plt.pie(counts)
            plt.legend(patches, labelnames, bbox_to_anchor=(0.1, 1), fontsize=15)
            plt.savefig(os.path.join(out_dir, "L1_" + self.tissue + "_quality_pie.png"))
