from typing import *
import os
import luigi
import development_mouse as dm


class WorkersTest(luigi.Task):

    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, f"number_workers.txt"))

    def run(self) -> None:
        with self.output().temporary_path() as out_file:
            f = open(out_file, "w")
            f.write(str(self.get_params()))
            f.close()
