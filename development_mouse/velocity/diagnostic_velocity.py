from typing import *
import os
import loompy
from scipy import sparse
import numpy as np
import cytograph as cg
import development_mouse as dm
import velocyto as vcy
import matplotlib.pyplot as plt
import luigi


class DiagnosticVelocity(luigi.Task):
    """Luigi Task to plot some diagnostic plots on the spliced, unspliced, ambigous molecules
    """
    tissue = luigi.Parameter(description="name of the tissue from tool specification file.",
                             always_in_help=True)

    def requires(self) -> List[luigi.Task]:
        """
        Arguments
        ---------
        `ClusterL1`:
            passing ``tissue``
        `AggregateL1`:
            passing ``tissue``
        `ExportL1`:
            passing ``tissue``
        """
        # NOTE: not sure it needs AggregateL1
        return [dm.ClusterL1(tissue=self.tissue),
                dm.AggregateL1(tissue=self.tissue),
                dm.ExportL1(tissue=self.tissue)]

    def output(self) -> luigi.Target:
        """
        Returns
        -------
        file: ``L1_[TISSUE]_exported"/velocity_{self.tissue}_SAU_fractions.pdf``
        """
        return luigi.LocalTarget(os.path.join(dm.paths().build, f"L1_{self.tissue}_exported", f"velocity_{self.tissue}_SAU_fractions.pdf"))

    def run(self) -> None:
        """Reads the output of `AggregateL1` and plots:
               - The fraction of splice, ambiguous, unspliced as ``velocity_[TISSUE]_exported/velocity_[TISSUE]_SAU_fractions.pdf``
        """
        logging = cg.logging(self, True)
        with self.output().temporary_path() as out_file:
            logging.info("Loading loom file in memory as a VelocytoLoom object")
            vlm = vcy.VelocytoLoom(self.input()[0].fn)
            logging.info("Plotting report on spliced, ambiguous, unpliced fraction")
            vlm.plot_fractions(save2file=out_file)
            # NOTE: other diagnostic plots if desired, but substitutet hte file for a path
