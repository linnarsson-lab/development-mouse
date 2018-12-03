from typing import *
import os
import csv
import pickle
import loompy
import numpy as np
import cytograph as cg
import development_mouse as dm
import luigi
import scipy.cluster.hierarchy as hierarchy
import numpy_groupies.aggregate_numpy as npg
import scipy.cluster.hierarchy as hc


class AggregateL1(luigi.Task):
	"""
	Aggregate all clusters in a new Loom file

	Parameters
	----------
	tissue: str
		name of the tissue from tool specification file
	n_markers: int, default=10
		the number of markers per cluster in the marker heatmap
	n_auto_genes: int, default=6
		number of genes to include in the auto-auto-annotation

	Raises
	------
	`ClusterL1`
		with the parameter tissue

	Returns
	-------
	On the output of `ClusterL1`:
	- Runs the `cytograph.Aggregator` and outputs thte results to itself (Note this breaks the luigi philosophy)
	- Runs the `cytograph.AutoAnnotator` and outputs the results to the ``L1_[TISSUE].agg.loom`` file
	- Runs the  `cytograph.AutoAutoAnnotator` and outputs the results to the ``L1_[TISSUE].agg.loom`` file

	Yields
	------
	file: ``L1_[TISSUE].agg.loom``

	"""
	tissue = luigi.Parameter()
	n_markers = luigi.IntParameter(default=10)
	n_auto_genes = luigi.IntParameter(default=6)

	def requires(self) -> List[luigi.Task]:
		return dm.ClusterL1(tissue=self.tissue)

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(dm.paths().build, f"L1_{self.tissue}.agg.loom"))

	def run(self) -> None:
		logging = cg.logging(self)
		with self.output().temporary_path() as out_file:
			ds = loompy.connect(self.input().fn)
			cg.Aggregator().aggregate(ds, out_file) # cg.Aggregator(self.n_markers).aggregate(ds, out_file) causes error
			dsagg = loompy.connect(out_file)

			logging.info("Computing auto-annotation")
			aa = cg.AutoAnnotator(root=dm.paths().autoannotation)
			aa.annotate_loom(dsagg)
			aa.save_in_loom(dsagg)

			logging.info("Computing auto-auto-annotation")
			n_clusters = dsagg.shape[1]
			(selected, selectivity, specificity, robustness) = cg.AutoAutoAnnotator(n_genes=self.n_auto_genes).fit(dsagg)
			dsagg.set_attr("MarkerGenes", np.array([" ".join(ds.row_attrs["Gene"][selected[:, ix]]) for ix in np.arange(n_clusters)]), axis=1)
			ds.close()
			np.set_printoptions(precision=1, suppress=True)
			dsagg.set_attr("MarkerSelectivity", np.array([str(selectivity[:, ix]) for ix in np.arange(n_clusters)]), axis=1)
			dsagg.set_attr("MarkerSpecificity", np.array([str(specificity[:, ix]) for ix in np.arange(n_clusters)]), axis=1)
			dsagg.set_attr("MarkerRobustness", np.array([str(robustness[:, ix]) for ix in np.arange(n_clusters)]), axis=1)
			dsagg.close()
