import os
from typing import *
import numpy as np
import scanpy.api as sc
import luigi
import loompy


class RunAGA(luigi.Task):
    """
    Luigi task to run scanpy AGA
    """
    tissue = luigi.Parameter(description="name of the tissue from tool specification file.",
                             always_in_help=True)

    def requires(self) -> List[luigi.Task]:
        return dm.ClusterL1(tissue=self.tissue)

    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, f"graph_{self.tissue}.h5"))

    def run(self) -> None:
        #add logging?
        with self.output().temporary_path()[:-3] as out_file:
            ds = loompy.connect(self.input())
            adata = sc.AnnData(np.transpose(ds[:, :]))
            adata.smp['Clusters'] = ds.col_attrs['Clusters'].astype('str')
            adata.var_names = ds.row_attrs['Gene']
            sc.pp.recipe_zheng17(adata, plot=True)  # replace with custom filtering?
            sc.tl.tsne(adata, n_jobs=16)
            ax = sc.pl.tsne(adata, color='Clusters')  # bug? why need to plot for tl.aga?
            sc.tl.aga(adata, node_groups='Clusters', n_jobs=16)
            sc.write(out_file, adata, ext='h5')