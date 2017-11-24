from typing import *
import luigi
import development_mouse as dm
import cytograph as cg


class AGAAnalysis(luigi.WrapperTask):
    """
    Luigi Wrapper Task that plots AGA for each tissue in the project using PlotAGA
    """
    project = luigi.Parameter(default="Development", description="str, default=Development\nSpecify set of ``tissues`` that will analyzed. Info is gathered using the ``pool_specification.tab`` file")

    def requires(self) -> Iterator[luigi.Task]:
        tissues = cg.PoolSpec().tissues_for_project(self.project)
        for tissue in tissues:
            yield dm.PlotAGA(tissue=tissue)
