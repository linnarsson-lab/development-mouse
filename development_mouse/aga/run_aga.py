import os
from typing import *
import numpy as np
import luigi
import loompy
import development_mouse as dm
import cytograph as cg
import logging


class RunAGA(luigi.Task):
    """
    Luigi task to run scanpy AGA
    """
    tissue = luigi.Parameter(description="name of the tissue from tool specification file.",
                             always_in_help=True)

    def requires(self) -> List[luigi.Task]:
        return dm.ClusterL1(tissue=self.tissue)

    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, f"AGA_{self.tissue}.h5"))

    def run(self) -> None:
        with self.output().temporary_path() as out_file:
            import scanpy.api as sc
            logging.info("Loading loom file into AnnData object")
            ds = loompy.connect(self.input().fn)
            adata = sc.AnnData(np.transpose(ds[:, :]))
            adata.smp['Clusters'] = ds.col_attrs['Clusters'].astype('str')
            adata.var_names = ds.row_attrs['Gene']
            logging.info("Filtering and t-SNE")
            sc.pp.recipe_zheng17(adata, plot=True)  # replace with custom filtering?
            sc.tl.tsne(adata, n_jobs=16)
            logging.info("Running AGA")
            sc.tl.aga(adata, groups='Clusters', n_jobs=16)
            sc.settings.writedir = ""
            logging.info("Writing AnnData to file")  
            sc.write(out_file, adata)
            os.rename(out_file + '.h5', out_file)
