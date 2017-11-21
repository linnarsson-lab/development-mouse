import numpy as np
import scanpy.api as sc
import luigi

class RunAGA(luigi.Task):
    """
    Luigi task to run scanpy's tl.aga 
    """
    def requires(self) -> :
        return dm.ClusterL1(tissue=self.tissue) 

    def output(self) -> :
        return luigi.LocalTarget(os.path.join(dm.paths().build, f"graph_{self.tissue}.h5"))

    def run(self) -> :
        with self.output().temporary_path() as out_file:

            sc.write(out_file, )
