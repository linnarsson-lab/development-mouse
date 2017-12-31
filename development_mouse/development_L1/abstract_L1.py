from typing import *
import os
import loompy
import numpy as np
import cytograph as cg
import development_mouse as dm
import luigi
import shutil


class AbstractL1(luigi.Task):
    """
    Perform graph abstraction on the clusters of ClusterL1

    Parameters
    ----------
    tissue: str
        name of the tissue from tool specification file
    confidence: float, default=0.5
        edges with confidence below `confidence` are removed

    Raises
    ------
    `ClusterL1` and `AggregateL1`
        with the parameter tissue

    Returns
    -------
    Run GraphAbstraction("simple").abstract on the clustered data

    Yields
    ------
    file: ``L1_[TISSUE]_abstract.agg.loom``

    """
    tissue = luigi.Parameter()
    confidence = luigi.FloatParameter(default=0.5)

    def requires(self) -> List[luigi.Task]:
        return [dm.ClusterL1(tissue=self.tissue), dm.AggregateL1(tissue=self.tissue)]

    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, f"L1_{self.tissue}_abstract.agg.loom"))

    def run(self) -> None:
        logging = cg.logging(self)
        with self.output().temporary_path() as out_file:
            ds = loompy.connect(self.input()[0].fn)
            shutil.copy(self.input()[1].fn, out_file)
            dsagg = loompy.connect(out_file)
            abstraction = cg.GraphAbstraction(kind="simple")
            agraph = abstraction.abstract(ds=ds, thresh_confid=self.confidence, unidirectional=False)
            ds.close()
            dsagg.set_edges("AbstractGraph", agraph.row, agraph.col, agraph.data)
            dsagg.close()
