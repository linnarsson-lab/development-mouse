from typing import *
import os
import csv
import numpy as np
import pandas as pd
import pickle
import logging
import luigi
import cytograph as cg
import development_mouse as dm
import loompy
import numpy.core.defchararray as npstr


def ixs_thatsort_a2b(a: np.ndarray, b: np.ndarray, check_content: bool=True) -> np.ndarray:
    "This is super duper magic sauce to make the order of one list to be like another"
    if check_content:
        assert len(np.intersect1d(a, b)) == len(a), f"The two arrays are not matching"
    return np.argsort(a)[np.argsort(np.argsort(b))]


class PoolL2(luigi.Task):
    """Luigi Task to pool the leaves from punchcards tree
    """
    punchcard_deck = dm.PunchcardParser()

    def requires(self) -> List[luigi.Task]:
        leaves: List[str] = self.punchcard_deck.prune_leaves()
        return [dm.Punchcard(card=l).requires() for l in leaves]

    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, "PoolL2.loom"))

    def run(self) -> None:
        # Pooling the files
        with self.output().temporary_path() as out_file:
            dsout: loompy.LoomConnection = None
            cluster_counter: int = 0
            reference_accession = None
            for punchcard in self.input():
                clusterP, exportP, *_ = punchcard  # It is confusing but it seems to not need .requires()
                ds = loompy.connect(clusterP.fn)
        
                reference_accession = ds.row_attrs["Accession"]
                order = ixs_thatsort_a2b(ds.row_attrs["Accession"], reference_accession)

                assert np.all(ds.col_attrs["Clusters"] != -1), "Some clusters are labeled -1 PoolL2 does not support that"
                # NOTE: This could be done in much easier way with loompy2
                for (ix, selection, vals) in ds.batch_scan_layers(axis=1, batch_size=dm.memory().axis1):
                    m = {}
                    for layer_name, chunk_of_matrix in vals.items():
                        m[layer_name] = vals[layer_name][order, :]  # NOTE: I think this is not needed since there is no selection  [:, selection - ix]
                    # NOTE: I don't need to do it this way I can do it outside the loop
                    ca = {}
                    for key in ds.col_attrs:
                        if key == "Clusters":
                            # NOTE Special attention not to merge clusters
                            ca["Clusters_original"] = ds.col_attrs[key][selection - ix]
                            ca["Clusters"] = ds.col_attrs[key][selection - ix] + cluster_counter
                        else:
                            ca[key] = ds.col_attrs[key][selection - ix]
                    ca["SourceFileName"] = np.full(len(selection), os.path.basename(clusterP.fn))

                    # Add data to the loom file
                    if dsout is None:
                        # create using main layer
                        dsout = loompy.create(out_file, m[""], ds.row_attrs, ca)
                        # Add layers
                        for layer_name, chunk_of_matrix in m.items():
                            if layer_name == "":
                                continue
                            dsout.set_layer(layer_name, chunk_of_matrix, dtype=chunk_of_matrix.dtype)
                    else:
                        dsout.add_columns(m, ca)
                ds.close()
                cluster_counter = np.max(dsout.col_attrs["Clusters"]) + 1
            dsout.close()
