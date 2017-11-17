from typing import *
import luigi
import multiprocessing


class threads(luigi.Config):
    limit = luigi.IntParameter(default=multiprocessing.cpu_count())
