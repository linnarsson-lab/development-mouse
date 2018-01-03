from typing import *
import os
import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import development_mouse as dm
import luigi


class Level1(luigi.WrapperTask):
    """Luigi Task to run a subset of level 1 Analysis
    """

    project = luigi.Parameter(default="Development")
    target = luigi.Parameter(default="All")  # one between Cortex, AllForebrain, ForebrainDorsal, ForebrainVentrolateral, ForebrainVentrothalamic, Midbrain, Hindbrain
    time = luigi.Parameter(default="E7-E18")  # later more specific autoannotation can be devised

    def requires(self) -> Iterator[luigi.Task]:
        tissues = dm.targets_map[self.target]
        for tissue in tissues:
            if dm.time_check(tissue, self.time):
                yield [dm.ClusterL1(tissue=tissue), dm.ExportL1(tissue=tissue)]  # NOTE this breaks all the the processes pipeline
