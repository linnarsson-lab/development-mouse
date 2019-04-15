from typing import *
import os
import logging
import loompy
import numpy as np
import cytograph as cg
import development_mouse as dm
import luigi


class AggregateL2(luigi.Task):  # Status: Ok
    """
    Summary statistics of all clusters in a new Loom file
    """
    n_markers = luigi.IntParameter(default=10, description="Number of markers considered by the Aggergator")
    n_auto_genes = luigi.IntParameter(default=6, description="Number of genes to use in the AutoAutoannotation")

    def requires(self) -> List[luigi.Task]:
        return dm.PoolL2()

    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, "PoolL2.agg.loom"))

    def run(self) -> None:
        logging = cg.logging(self)
        with self.output().temporary_path() as out_file:
            ds = loompy.connect(self.input().fn)
            logging.info("Aggregating L2 Pool")
            agg_spec = {
                        "Age": "tally",
                        "Clusters": "first",
                        "Class": "mode",
                        "_Total": "mean",
                        "Sex": "tally",
                        "Tissue": "tally",
                        "SampleID": "tally",
                        "TissuePool": "first",
                        "Outliers": "mean",
                        "Clusters_original": "mode",
                        "SourceFileName": "mode"
			            }
            cg.Aggregator().aggregate(ds, out_file, agg_spec=agg_spec)
            # cg.Aggregator(self.n_markers).aggregate(ds, out_file, batch_size=dm.memory().axis0)
            dsagg = loompy.connect(out_file)

            logging.info("Computing auto-annotation of L2 Pool Aggregated")
            aa = cg.AutoAnnotator(root=dm.paths().autoannotation)
            aa.annotate_loom(dsagg)
            aa.save_in_loom(dsagg)
            
            ds.close()
            # logging.info("Computing auto-auto-annotation")
            # n_clusters = dsagg.shape[1]
            # (selected, selectivity, specificity, robustness) = cg.AutoAutoAnnotator(n_genes=self.n_auto_genes).fit(dsagg)
            # dsagg.set_attr("MarkerGenes", np.array([" ".join(ds.row_attrs["Gene"][selected[:, ix]]) for ix in np.arange(n_clusters)]), axis=1)
            # ds.close()
            # np.set_printoptions(precision=1, suppress=True)
            # dsagg.set_attr("MarkerSelectivity", np.array([str(selectivity[:, ix]) for ix in np.arange(n_clusters)]), axis=1)
            # dsagg.set_attr("MarkerSpecificity", np.array([str(specificity[:, ix]) for ix in np.arange(n_clusters)]), axis=1)
            # dsagg.set_attr("MarkerRobustness", np.array([str(robustness[:, ix]) for ix in np.arange(n_clusters)]), axis=1)
            # dsagg.close()
