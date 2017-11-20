from typing import *
import os
import csv
import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import development_mouse as dm
import luigi
from collections import defaultdict


class Punchcard(luigi.WrapperTask):  # Status: check what it should return
	"""
	Luigi Task Wrapper to run a punchcard

	`card` needs to match th name specified in the .yaml file in the folder ../punchcards
	"""
	
	card = luigi.Parameter(description="Name of the punchcard to run")

	def requires(self) -> List[List[luigi.Task]]:
		punchcard_obj = dm.PunchcardParser()[self.card]
		# To make it more general we can avoid to have ExportPunchcard as the first instance
		other_tasks = []
		for task in cg.parse_analysis_todo(punchcard_obj):
			other_tasks.append(task(card=self.card))
		return [dm.ExportAnalysis(card=self.card), *other_tasks]  # Not sure why but before t was [[]]
