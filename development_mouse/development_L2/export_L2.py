from typing import *
import os
import logging
import loompy
from scipy import sparse
import numpy as np
import cytograph as cg
import development_mouse as dm
import luigi


class ExportL2(luigi.Task):
    """
    Luigi Task to export summary files for a punchcard analysis
    """
    n_markers = luigi.IntParameter(default=10, description="number of markers to export")

    def requires(self) -> luigi.Task:
        # NOTE before the order was AggregatePunchcard, ClusterPunchcard but did not make sense
        return {"PoolL2": dm.PoolL2(),
                "AggregateL2": dm.AggregateL2()}

    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, f"PoolL2_exported"))

    def run(self) -> None:
        logging = cg.logging(self, True)
        logging.info("Exporting cluster data")
        with self.output().temporary_path() as out_dir:
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            ds = loompy.connect(self.input()["PoolL2"].fn)
            dsagg = loompy.connect(self.input()["AggregateL2"].fn)
            
            # dsagg.export(os.path.join(out_dir, "PoolL2_expression.tab"))
            # dsagg.export(os.path.join(out_dir, "PoolL2_enrichment.tab"), layer="enrichment")
            # dsagg.export(os.path.join(out_dir, "PoolL2_enrichment_q.tab"), layer="enrichment_q")
            # dsagg.export(os.path.join(out_dir, "PoolL2_trinaries.tab"), layer="trinaries")

            # dm.plot_dendrogram(dsagg)
            # ds = loompy.connect(self.requires().input().fn)
            # logging.info(f"Plotting marker heatmap for {self.card}")
            # cg.plot_markerheatmap(ds, dsagg, n_markers_per_cluster=self.n_markers, out_file=os.path.join(out_dir, f"{self.card}_heatmap.pdf"))

            logging.info("Plotting manifold graph with auto-annotation")
            tags = list(dsagg.col_attrs["AutoAnnotation"])
            cg.plot_graph(ds, os.path.join(out_dir, "L2_manifold.aa.png"), tags)

            logging.info("Plotting manifold graph with auto-auto-annotation")
            tags = list(dsagg.col_attrs["MarkerGenes"])
            cg.plot_graph(ds, os.path.join(out_dir, "L2_manifold.aaa.png"), tags)

            logging.info("Plotting marker heatmap")
            cg.plot_markerheatmap(ds, dsagg, n_markers_per_cluster=self.n_markers, out_file=os.path.join(out_dir, "L2_heatmap.pdf"))

            logging.info("Plotting UMAP")
            cg.plot_graph(ds, os.path.join(out_dir, "L2_UMAP_manifold.aaa.png"), tags, embedding="UMAP")
            logging.info("Plotting UMI and gene counts")
            cg.plot_umi_genes(ds, out_file=os.path.join(out_dir, "L2_umi_genes.png"))
            logging.info("Plotting factors")
            cg.plot_factors(ds, base_name=os.path.join(out_dir, "L2_factors"))
            logging.info("Plotting cell cycle")
            cg.plot_cellcycle(ds, out_file=os.path.join(out_dir, "L2_cellcycle.png"))
            logging.info("Plotting markers")
            cg.plot_markers(ds, out_file=os.path.join(out_dir, "L2_markers.png"))
            logging.info("Plotting neighborhood diagnostics")
            cg.plot_radius_characteristics(ds, out_file=os.path.join(out_dir, "L2_neighborhoods.png"))
            logging.info("Plotting batch covariates")
            cg.plot_batch_covariates(ds, out_file=os.path.join(out_dir, "L2_batches.png"))
            cg.ClusterValidator().fit(ds, os.path.join(out_dir, "L2_cluster_pp.png"))
            logging.info("Plotting embedded velocity")
            cg.plot_embedded_velocity(ds, out_file=os.path.join(out_dir, "L2_velocity.png"))
            logging.info("Plotting TFs")
            cg.plot_TFs(ds, dsagg, out_file_root=os.path.join(out_dir, "L2"))