from typing import *
import os
import scanpy.api as sc
import luigi
import development_mouse as dm
import matplotlib.pyplot as plt


class PlotAGA(luigi.Task):
    """
    Luigi task to plot scanpy AGA
    """
    tissue = luigi.Parameter(description="name of the tissue from tool specification file.",
                             always_in_help=True)

    def requires(self) -> List[luigi.Task]:
        return dm.RunAGA(tissue=self.tissue)

    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, f"AGA_{self.tissue}_graph.png"))

    def run(self) -> None:
        with self.output().temporary_path() as out_file:
            adata = sc.read(self.input().fn)
            print('outfile:', self.input().fn)            
            sc.pl.aga(adata, color='Clusters', layout='fr')
            plt.savefig(out_file)
