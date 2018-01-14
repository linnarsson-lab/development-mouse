from typing import *
import os
from shutil import copyfile
import numpy as np
import logging
import luigi
import cytograph as cg
import development_mouse as dm
import loompy
import logging
from scipy import sparse
from sklearn.neighbors import KNeighborsClassifier
import pickle
import networkx as nx
import hdbscan
import gc
from sklearn.cluster import DBSCAN


class ClusterPunchcard(luigi.Task):  # Status: OK
    """
    Clustering for a Punchcard pool
    """
    card = luigi.Parameter(description="Name of the punchcard to use")
    n_genes = luigi.IntParameter(default=1000, description="""(default=1000) The number of genes used in manifold learning""")
    manifold_learning = luigi.IntParameter(default=1, description="int, default=1\nWhether to use `cytograph.ManifoldLearning` or `cytograph.ManifoldLearning2`", always_in_help=True)
    gtsne = luigi.BoolParameter(default=True, description="(default=True) Use graph t-SNE for layout")
    alpha = luigi.FloatParameter(default=1, description="(default=1) The scale parameter for multiscale KNN")
    filter_geneset = luigi.Parameter(default=None, description="The path of a file containing as rows the gene symbol of genes to excluded (Note: despite the name it can be used for any gene set)")
    layer = luigi.Parameter(default=None, description="Layer used for manifold learning (i.e. the matrix used to compute PCA).Currently it only has effects when using `cytograph.ManifoldLearning` and not `cytograph.ManifoldLearning2`")

    def requires(self) -> luigi.Task:
        return dm.PunchcardPool(card=self.card)

    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, f"{self.card}.loom"))

    def run(self) -> None:
        logging = cg.logging(self)
        with self.output().temporary_path() as out_file:
            ds = loompy.connect(self.input().fn)
            dsout: loompy.LoomConnection = None

            logging.info("Removing invalid cells")
            for (ix, selection, vals) in ds.batch_scan_layers(cells=np.where(ds.col_attrs["_Valid"] == 1)[0], layers=ds.layer.keys(), batch_size=dm.memory().axis1, axis=1):
                ca = {key: val[selection] for key, val in ds.col_attrs.items()}
                if dsout is None:
                    # NOTE Loompy Create should support multilayer !!!!
                    if type(vals) is dict:
                        dsout = loompy.create(out_file, vals[""], row_attrs=ds.row_attrs, col_attrs=ca, dtype=vals[""].dtype)
                        for layername, layervalues in vals.items():
                            if layername != "":
                                dsout.set_layer(layername, layervalues, dtype=layervalues.dtype)
                        dsout = loompy.connect(out_file)
                    else:
                        loompy.create(out_file, vals, row_attrs=ds.row_attrs, col_attrs=ca)
                        dsout = loompy.connect(out_file)
                else:
                    dsout.add_columns(vals, ca)
            # dsout.close() causes an exception; disabling gc fixes it. See https://github.com/h5py/h5py/issues/888
            gc.disable()
            dsout.close()
            gc.enable()

            if self.manifold_learning == 2:
                logging.info("Learning the manifold")
                ds = loompy.connect(out_file)
                ml = cg.ManifoldLearning2(n_genes=self.n_genes, gtsne=self.gtsne, alpha=self.alpha, filter_cellcycle=self.filter_geneset, layer=self.layer)
                (knn, mknn, tsne) = ml.fit(ds)
                ds.set_edges("KNN", knn.row, knn.col, knn.data, axis=1)
                ds.set_edges("MKNN", mknn.row, mknn.col, mknn.data, axis=1)
                ds.set_attr("_X", tsne[:, 0], axis=1)
                ds.set_attr("_Y", tsne[:, 1], axis=1)

                logging.info("Clustering on the manifold")
                cls = cg.Clustering(method="mknn_louvain", min_pts=10)
                labels = cls.fit_predict(ds)
                ds.set_attr("Clusters", labels, axis=1)
                logging.info(f"Found {labels.max() + 1} clusters")
                cg.Merger(min_distance=0.2).merge(ds)
                logging.info(f"Merged to {ds.col_attrs['Clusters'].max() + 1} clusters")
                ds.close()
            elif self.manifold_learning == 3:
                logging.info("Learning the manifold")
                ds = loompy.connect(out_file)
                ml = cg.ManifoldLearning2(n_genes=self.n_genes, gtsne=self.gtsne, alpha=self.alpha, filter_cellcycle=self.filter_geneset, layer=self.layer)
                (knn, mknn, tsne) = ml.fit(ds)
                ds.set_edges("KNN", knn.row, knn.col, knn.data, axis=1)
                ds.set_edges("MKNN", mknn.row, mknn.col, mknn.data, axis=1)
                ds.set_attr("_X", tsne[:, 0], axis=1)
                ds.set_attr("_Y", tsne[:, 1], axis=1)

                logging.info("Clustering on the manifold")
                pl = cg.PolishedLouvain()
                # This is for loompy1 For loompy2 just: labels = pl.fit_predict(dsout.col_graphs.MKNN, tsne)
                a, b, w = ds.get_edges("MKNN", axis=1)
                knn = sparse.coo_matrix((w, (a, b)), shape=(ds.shape[1], ds.shape[1]))
                labels = pl.fit_predict(knn, tsne)

                ds.set_attr("Clusters", labels + 1, axis=1)
                ds.set_attr("Outliers", (labels == -1).astype('int'), axis=1)
                logging.info(f"Found {labels.max() + 1} clusters")
                dsout.close()
            else:
                dsout = loompy.connect(out_file)
                ml = cg.ManifoldLearning(n_genes=self.n_genes, gtsne=self.gtsne, alpha=self.alpha, filter_cellcycle=self.filter_geneset, layer=self.layer)
                (knn, mknn, tsne) = ml.fit(dsout)

                dsout.set_edges("KNN", knn.row, knn.col, knn.data, axis=1)
                dsout.set_edges("MKNN", mknn.row, mknn.col, mknn.data, axis=1)
                dsout.set_attr("_X", tsne[:, 0], axis=1)
                dsout.set_attr("_Y", tsne[:, 1], axis=1)

                min_pts = 10
                eps_pct = 90
                cls = cg.Clustering(method=dm.cluster().method, outliers=not dm.cluster().no_outliers, min_pts=min_pts, eps_pct=eps_pct)
                labels = cls.fit_predict(dsout)
                dsout.set_attr("Clusters", labels, axis=1)
                n_labels = np.max(labels) + 1
                dsout.close()
