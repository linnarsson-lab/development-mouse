from typing import *
import os
from shutil import copyfile
import numpy as np
import luigi
import gc
import cytograph as cg
import development_mouse as dm
import loompy
import logging
import pickle
from scipy import sparse
from scipy.special import polygamma
from sklearn.cluster import AgglomerativeClustering, KMeans, Birch
from sklearn.decomposition import PCA, IncrementalPCA, FastICA
from sklearn.manifold import TSNE
from sklearn.metrics import pairwise_distances
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.neighbors import BallTree, NearestNeighbors, kneighbors_graph, KNeighborsClassifier
from sklearn.preprocessing import scale
from sklearn.svm import SVR
from scipy.stats import ks_2samp
import networkx as nx
import hdbscan
from sklearn.cluster import DBSCAN


class ClusterL1(luigi.Task):
    """Level 1 Clustering
    """

    tissue = luigi.Parameter(description="Name of the tissue from tool specification file")
    n_genes = luigi.IntParameter(default=1000, description="""(default=1000) The number of genes used in manifold learning""")
    manifold_learning = luigi.IntParameter(default=1, description="(default=1) Whether to use `cytograph.ManifoldLearning` or `cytograph.ManifoldLearning2`")
    gtsne = luigi.BoolParameter(default=True, description="(default=True) Use graph t-SNE for layout")
    alpha = luigi.FloatParameter(default=1, description="(default=1) The scale parameter for multiscale KNN")
    filter_geneset = luigi.Parameter(default="None", description="(path) The path of a file containing as rows the gene symbol of genes to excluded (Note: despite the name it can be used for any gene set)")
    layer = luigi.Parameter(default="None", description="Layer used for manifold learning (i.e. the matrix used to compute PCA). Currently it only has effects when using `cytograph.ManifoldLearning` and not `cytograph.ManifoldLearning2`")

    def requires(self) -> luigi.Task:
        return {"PrepareTissuePool": dm.PrepareTissuePool(tissue=self.tissue)}

    def output(self) -> luigi.Target:
        """
        Returns
        -------
        file: ``L1_[TISSUE].loom``

        Reads the output of `PrepareTissuePool` and does:
            - runs the `cytograph.batch_scan_layers` to remove cells that are not valid (almost creating a copy)
            - runs `cytograph.ManifoldLearning` or `cytograph.ManifoldLearning2`
            - runs `cytograph.Clustering` on the learned mainfild
            - if `cytograph.ManifoldLearning2` also runs `cytograph.Merger`
        """
        return luigi.LocalTarget(os.path.join(dm.paths().build, "L1_" + self.tissue + ".loom"))

    def run(self) -> None:
        """
        Other Parameters
        ----------------
        cluster_method: `cluster_config.cluster.method`, default="dbscan"
            Select the method for clustering. Valid: "hdbscan", "dbscan", "multilev", "wmultilev", "mknn_louvain", "louvain" (same as None)
        cluster_no_outliers: `cluster_config.cluster.no_outliers`, default=False
            Whether to consider in the clustering cells that have been marked as outliers in `PrepareTissuePool`
        memory_batch: `memory_config.memory.batch`, default=2000
            The size of the batches that are used by `cytograph.batch_scan_layers`
        """
        if self.filter_geneset == "None":
            self.filter_geneset = None
        if self.layer == "None":
            self.layer = None
        logging = cg.logging(self)
        with self.output().temporary_path() as out_file:
            ds = loompy.connect(self.input()["PrepareTissuePool"].fn)
            dsout: loompy.LoomConnection = None

            #logging.info("Deserializing QC Classifier")
            #knc: KNeighborsClassifier = pickle.load(open(os.path.join(self.input()["MakeQualityClassifier"].fn, "QC_Classifier.pickle"), "rb"))
            #logging.info("Reading NameQualityCluster file")
            # cluster_mapping = {int(i.split(":")[0]): i.split(":")[1] for i in open(self.input()["NameQualityClusters"].fn).read().rstrip().split()}
            #initial_cell_size = ds.col_attrs["SplicedTotal"]
            #initial_Ucell_size = ds.col_attrs["UnsplicedTotal"]
            #detected_genes = ds.col_attrs["TotalMolNoAmbiguous"]
            #mito_size = ds.col_attrs["MitocondrialTotal"]
            #ribo_size = ds.col_attrs["RibosomalTotal"]
            #X = np.column_stack((initial_cell_size, initial_Ucell_size, detected_genes, mito_size, ribo_size))
            #X_log = np.log2(X + 1)

            #logging.info("Using the QC Classifier to set QualityClass")
            #predicted = knc.predict(X_log)
            # qc_luster_labels = np.array([cluster_mapping[i] for i in predicted])
            #ds.set_attr(name="QualityClass", values=predicted, axis=1)

            # NOTE for now the quality class is only written and not used anywhere

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
                # Note dm.cluster().no_outliers is used only in the clustering dbscan and mknn_louvain
                cls = cg.Clustering(method=dm.cluster().method, outliers=not dm.cluster().no_outliers, min_pts=min_pts, eps_pct=eps_pct)
                labels = cls.fit_predict(dsout)
                dsout.set_attr("Clusters", labels, axis=1)
                if np.any(labels == -1):
                    ds.set_attr("Outliers", (labels == -1).astype('int'), axis=1)
                # n_labels = np.max(labels) + 1
                dsout.close()
