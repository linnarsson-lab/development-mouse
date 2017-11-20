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


class ExportPunchcard(luigi.Task):
    """
    Luigi Task to export summary files for a punchcard analysis
    """
    card = luigi.Parameter(description="Name of punchcard")
    n_markers = luigi.IntParameter(default=10, description="number of markers to export")

    def requires(self) -> List[luigi.Task]:
        return [dm.AggregatePunchcard(analysis=self.card),
                dm.ClusterPunchcard(analysis=self.card)]

    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, f"{self.card}_exported"))

    def run(self) -> None:
        logging = cg.logging(self, True)
        logging.info("Exporting cluster data")
        with self.output().temporary_path() as out_dir:
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            dsagg = loompy.connect(self.input()[0].fn)
            dsagg.export(os.path.join(out_dir, f"{self.card}_expression.tab"))
            dsagg.export(os.path.join(out_dir, f"{self.card}_enrichment.tab"), layer="enrichment")
            dsagg.export(os.path.join(out_dir, f"{self.card}_enrichment_q.tab"), layer="enrichment_q")
            dsagg.export(os.path.join(out_dir, f"{self.card}_trinaries.tab"), layer="trinaries")

            ds = loompy.connect(self.input()[1].fn)

            logging.info(f"Plotting manifold graph with auto-annotation for {self.card}")
            tags = list(dsagg.col_attrs["AutoAnnotation"])
            cg.plot_graph(ds, os.path.join(out_dir, f"{self.card}_manifold.aa.png"), tags)

            logging.info(f"Plotting manifold graph with auto-auto-annotation for {self.card}")
            tags = list(dsagg.col_attrs["MarkerGenes"])
            cg.plot_graph(ds, os.path.join(out_dir, f"{self.card}_manifold.aaa.png"), tags)

            logging.info(f"Plotting marker heatmap for {self.card}")
            cg.plot_markerheatmap(ds, dsagg, n_markers_per_cluster=self.n_markers, out_file=os.path.join(out_dir, f"{self.card}_heatmap.pdf"))
