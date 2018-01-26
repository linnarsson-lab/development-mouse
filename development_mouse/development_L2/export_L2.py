from typing import *
import os
import logging
import loompy
from scipy import sparse
import numpy as np
import networkx as nx
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
        return dm.AggregateL2()

    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, f"PoolL2_exported"))

    def run(self) -> None:
        logging = cg.logging(self, True)
        logging.info("Exporting cluster data")
        with self.output().temporary_path() as out_dir:
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            dsagg = loompy.connect(self.input().fn)
            dsagg.export(os.path.join(out_dir, "PoolL2_expression.tab"))
            dsagg.export(os.path.join(out_dir, "PoolL2_enrichment.tab"), layer="enrichment")
            dsagg.export(os.path.join(out_dir, "PoolL2_enrichment_q.tab"), layer="enrichment_q")
            dsagg.export(os.path.join(out_dir, "PoolL2_trinaries.tab"), layer="trinaries")

            # ds = loompy.connect(self.requires().input().fn)
            # logging.info(f"Plotting marker heatmap for {self.card}")
            # cg.plot_markerheatmap(ds, dsagg, n_markers_per_cluster=self.n_markers, out_file=os.path.join(out_dir, f"{self.card}_heatmap.pdf"))
