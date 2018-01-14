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
    tissue = luigi.Parameter(description="name of the tissue from tool specification file.")

    n_markers = luigi.IntParameter(default=10,
                                   description="number of markers visualized in the heatmap")

    def requires(self) -> Dict[str, luigi.Task]:
        """
        Arguments
        ---------
        `AbstractL1`:
            passing ``tissue``
        `ClusterL1`:
            passing ``tissue``
        """
        return {f"AbstractL1(tissue={self.tissue})": dm.AbstractL1(tissue=self.tissue),
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
        Reads the output of `AbstractL1` and does:
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
            dsagg = loompy.connect(self.input()[f"AbstractL1(tissue={self.tissue})"].fn)

            dsagg.export(os.path.join(out_dir, f"L1_{self.tissue}_expression.tab"))
            dsagg.export(os.path.join(out_dir, f"L1_{self.tissue}_enrichment.tab"), layer="enrichment")
            dsagg.export(os.path.join(out_dir, f"L1_{self.tissue}_enrichment_q.tab"), layer="enrichment_q")
            dsagg.export(os.path.join(out_dir, f"L1_{self.tissue}_trinaries.tab"), layer="trinaries")
            # NOTE: maybe we should have dsagg.export(*, layer="abstraction")

            ds = loompy.connect(self.input()[f"ClusterL1(tissue={self.tissue})"].fn)

            logging.info("Plotting manifold graph with auto-annotation")
            tags = list(dsagg.col_attrs["AutoAnnotation"])
            cg.plot_graph(ds, os.path.join(out_dir, f"L1_{self.tissue}_manifold.aa.png"), tags)
            logging.info("Plotting manifold graph with auto-annotation, colored by age")
            cg.plot_graph_age(ds, os.path.join(out_dir, f"L1_{self.tissue}_manifold.age.png"), tags)

            logging.info("Plotting abstracted graph with auto-annotation")
            dm.plot_abs_graph(ds, dsagg, os.path.join(out_dir, f"L1_{self.tissue}_absgraph.aa.png"), tags)

            logging.info("Plotting manifold graph with auto-auto-annotation")
            tags = list(dsagg.col_attrs["MarkerGenes"])
            cg.plot_graph(ds, os.path.join(out_dir, f"L1_{self.tissue}_manifold.aaa.png"), tags)

            logging.info("Plotting marker heatmap")
            cg.plot_markerheatmap(ds, dsagg, n_markers_per_cluster=self.n_markers, out_file=os.path.join(out_dir, f"L1_{self.tissue}_heatmap.pdf"))

            logging.info("Plotting quality class on t-SNE")
            tags = list(dsagg.col_attrs["AutoAnnotation"])
            cluster_mapping = {int(i.split(":")[0]): i.split(":")[1] for i in open(self.input()["NameQualityClusters"].fn).read().rstrip().split()}
            dm.plot_quality_graph(ds, dsagg, out_file=os.path.join(out_dir, f"L1_{self.tissue}_quality_tsne.png"),
                                  cluster_mapping=cluster_mapping, tags=tags)
                    
            logging.info("Plotting quality class in pie chart")
            plt.figure(None, (10, 10))
            labels = ds.col_attrs["QualityClass"].astype(int)
            unique, counts = np.unique(labels, return_counts=True)
            labelnames = [cluster_mapping[ix] for ix in unique]
            patches, texts = plt.pie(counts)
            plt.legend(patches, labelnames, bbox_to_anchor=(0.1, 1), fontsize=15)
            plt.savefig(os.path.join(out_dir, "L1_" + self.tissue + "_quality_pie.png"))
