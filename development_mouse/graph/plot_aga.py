import numpy as np
import scanpy.api as sc
import luigi
import development_mouse as dm


class PlotAGA(luigi.Task):
    """
    Luigi task to plot scanpy AGA
    """
    tissue = luigi.Parameter(description="name of the tissue from tool specification file.",
                             always_in_help=True)

    def requires(self) -> List[luigi.Task]:
        return dm.RunAGA(tissue=self.tissue)

    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, f"AGA_{self.tissue}_export"))

    def run(self) -> :