from typing import *
import os
import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi


class Level1Analysis(luigi.WrapperTask):
	"""
	Luigi Wrapper Task to trigger all Level 1 analyses through ExportL1

	Arguments
	---------
	project: str, default=Development
		Specify set of ``tissues`` that will analyzed. Info is gathered using the ``pool_specification.tab`` file
	
	Returns
	-------
	Yields ``ExportL1(tissue) for tissue in project``
	"""

	project = luigi.Parameter(default="Development")  # NOTE: This could be hardcoded since the pipeline split, but I leave it as an option for more flexibility. 

	def requires(self) -> Iterator[luigi.Task]:
		tissues = cg.PoolSpec().tissues_for_project(self.project)
		for tissue in tissues:
			yield cg.ExportL1(tissue=tissue)
