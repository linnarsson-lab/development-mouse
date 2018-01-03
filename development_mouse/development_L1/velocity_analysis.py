from typing import *
import os
import logging
import cytograph as cg
import development_mouse as dm
import luigi


class VelocityAnalysis(luigi.WrapperTask):
    """
    Luigi Wrapper Task to trigger all Velocity analyses through VisualizeVelocity

    Parameters
    ----------
    project: str, default=Development
        Specify set of ``tissues`` that will analyzed. Info is gathered using the ``pool_specification.tab`` file
    
    Raises
    ------
    Yields ``VisualizeVelocity(tissue) for tissue in project``
    """

    project = luigi.Parameter(default="Development", description="str, default=Development\nSpecify set of ``tissues`` that will analyzed. Info is gathered using the ``pool_specification.tab`` file")
    # NOTE: This could be hardcoded since the pipeline split, but I leave it as an option for more flexibility.

    def requires(self) -> Iterator[luigi.Task]:
        tissues = cg.PoolSpec().tissues_for_project(self.project)
        for tissue in tissues:
            yield dm.VisualizeVelocity(tissue=tissue)
