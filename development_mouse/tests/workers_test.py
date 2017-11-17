from typing import *
import os
import luigi
import development_mouse as dm


class WorkersTest(luigi.Task):
    """Luigi Task to plot some diagnostic plots on the spliced, unspliced, ambigous molecules
    """

    def output(self) -> luigi.Target:
        """
        Returns
        -------
        file: ``L1_[TISSUE]_exported"/velocity_{self.tissue}_SAU_fractions.pdf``
        """
        return luigi.LocalTarget(os.path.join(dm.paths().build, f"number_workers.txt"))

    def run(self) -> None:
        """Reads the output of `AggregateL1` and plots:
               - The fraction of splice, ambiguous, unspliced as ``velocity_[TISSUE]_exported/velocity_[TISSUE]_SAU_fractions.pdf``
        """
        with self.output().temporary_path() as out_file:
            f = open(out_file, "w")
            f.write(str(self.workers))
            f.write(str(self.deps))
            f.close()
