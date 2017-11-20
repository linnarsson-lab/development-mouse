from typing import *
import os
import logging
import pickle
import loompy
import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import KNeighborsClassifier
import matplotlib as mpl
import matplotlib.pyplot as plt
import velocyto as vcy
import cytograph as cg
import development_mouse as dm
import luigi


class MakeQualityClassifier(luigi.Task):
    """
    Luigi Task to train a quality classifier using summary statistics such as number of spliced, unspliced, mitocondrial, ribosomal
    """

    project = luigi.Parameter(default="Development")
    components = luigi.IntParameter(default=15, description="Number of components of the ")
    n_neighbors = luigi.IntParameter(default=100, description="The number of neighbours to use to classify the quality of a cell")

    def requires(self) -> Iterator[luigi.Task]:
        tissues = cg.PoolSpec().tissues_for_project(self.project)
        for tissue in tissues:
            yield dm.PrepareTissuePool(tissue=tissue)

    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, "QC_classifier"))

    def run(self) -> None:
        with self.output().temporary_path() as path_name:
            if not os.path.exists(path_name):
                os.mkdir(path_name)
            logging.info("Get features to compute quality classifier")

            for input_obj in self.input():
                logging.info(f"Get summary stats from {input_obj.fn}")
                ds = loompy.connect(input_obj.fn)

                initial_cell_size = ds.col_attrs["SplicedTotal"]
                initial_Ucell_size = ds.col_attrs["UnsplicedTotal"]
                detected_genes = ds.col_attrs["TotalMolNoAmbiguous"]
                mito_size = ds.col_attrs["MitocondrialTotal"]
                ribo_size = ds.col_attrs["RibosomalTotal"]

                try:
                    sum_stats = np.column_stack((initial_cell_size, initial_Ucell_size, detected_genes, mito_size, ribo_size))
                    X = np.row_stack((X, sum_stats))
                except NameError:
                    X = np.column_stack((initial_cell_size, initial_Ucell_size, detected_genes, mito_size, ribo_size))
            logging.info("Summary stats loaded")
            feats_name = ["SplicedTotal", "UnsplicedTotal", "TotalMolNoAmbiguous", "MitocondrialTotal", "RibosomalTotal"]
            X_log = np.log2(X + 1)

            bics = []
            bics_avg = []
            bics_std = []
            gmms = []
            N = self.components
            for n_components in range(1, N):
                logging.info(f"Fitting a GMM with n_components={n_components} to the summary statistics")
                gmm_tmp = GaussianMixture(n_components=n_components, n_init=10, tol=5e-4, max_iter=150, random_state=19900715)
                ix = np.random.choice(X.shape[0], size=int(X.shape[0] * 1), replace=False)
                gmm_tmp.fit(X_log[ix, :])
                bics.append(gmm_tmp.bic(X_log[ix, :]))
                tmp_bics = []
                for rep in range(20):
                    ix = np.random.choice(X.shape[0], size=int(X.shape[0] * 0.25), replace=True)
                    tmp_bics.append(gmm_tmp.bic(X_log[ix, :]))
                bics_avg.append(np.mean(tmp_bics))
                bics_std.append(np.std(tmp_bics))
                gmms.append(gmm_tmp)

            bics_avg_a = np.array(bics_avg)
            bics_std_a = np.array(bics_std)
            for i in range(len(bics_avg)):
                if np.alltrue(bics_avg_a[i] + 1.96 * bics_std_a[i] < (bics_avg_a - 1.96 * bics_std_a)[0:i]):
                    current_min = i

            logging.info(f"Chosing n_components by BIC yields {np.argmin(bics)} components")
            logging.info(f"Chosing n_components by BIC bootstrapping yields {current_min} components")

            logging.info(f"Plot the BIC values for different number of components")
            plt.figure(None, (10, 10))
            plt.bar(np.arange(1, N), bics_avg_a, width=.2, yerr=bics_std_a)
            plt.bar(np.arange(1, N)[current_min:current_min + 1], bics_avg[current_min], width=.2, yerr=bics_std_a[current_min], color="r")
            plt.ylim((np.min(bics_avg_a) * 0.95, np.max(bics_avg_a) * 1.05))
            plt.savefig(os.path.join(path_name, "BIC_GMM_pick.png"))

            chosen_N = np.argmin(bics)
            gmm = gmms[chosen_N]  # np.argmin(bics) or  current_min
            labels = gmm.predict(X_log)

            colorandum = plt.cm.tab20(labels / np.max(labels))
            
            # subsample for plotting
            ix = np.random.choice(X.shape[0], size=min(100000, X.shape[0]), replace=False)

            logging.info("Ploting QC_Classes_GMM.png")
            plt.figure(None, (15, 15))
            gs = plt.GridSpec(X.shape[1], X.shape[1])
            for i in range(X.shape[1]):
                for j in range(X.shape[1]):
                    if i <= j:
                        pass
                    else:
                        ax = plt.subplot(gs[i, j])
                        sc = vcy.scatter_viz(X_log[ix, j], X_log[ix, i],
                                             s=3, alpha=0.06, c=labels[ix].astype(float),
                                             vmin=0, vmax=20, cmap=plt.cm.tab20)
                        plt.xlabel(feats_name[j])
                        plt.ylabel(feats_name[i])
                        for c in range(np.max(labels) + 1):
                            plt.text(np.median(X_log[ix, j][labels[ix] == c]), np.median(X_log[ix, i][labels[ix] == c]), f"{c}",
                                     color="k", fontdict={"size": 10},
                                     bbox=dict(facecolor='w', edgecolor='none', alpha=0.2))
            plt.tight_layout()
            plt.savefig(os.path.join(path_name, "QC_Classes_GMM.png"), dpi=300)

            plt.figure(None, (15, 15))
            c = 0
            for c in range(chosen_N):
                logging.info(f"Ploting QC_Class{c}_GMM.png")
                gs = plt.GridSpec(X.shape[1], X.shape[1])
                for i in range(X.shape[1]):
                    for j in range(X.shape[1]):
                        if i <= j:
                            pass
                        else:
                            ax = plt.subplot(gs[i, j])
                            vcy.scatter_viz(X_log[ix, j], X_log[ix, i], s=1, alpha=0.04, c="0.6")
                            sc = vcy.scatter_viz(X_log[ix, j][labels[ix] == c], X_log[ix, i][labels[ix] == c],
                                                 s=3, alpha=0.08, c=labels[ix][labels[ix] == c].astype(float), vmin=0, vmax=20, cmap=plt.cm.tab20)
                            plt.xlabel(feats_name[j])
                            plt.ylabel(feats_name[i])
                            plt.text(np.median(X_log[ix, j][labels[ix] == c]), np.median(X_log[ix, i][labels[ix] == c]), f"{c}",
                                     color="k", fontdict={"size": 10},
                                     bbox=dict(facecolor='w', edgecolor='none', alpha=0.2))
                                
                plt.tight_layout()
                plt.savefig(os.path.join(path_name, f"QC_Class{c}_GMM.png"), dpi=300)

            logging.info("Preparing (e.g. is training-less) the K Nearest Neighbors Classifier")
            knc = KNeighborsClassifier(n_neighbors=self.n_neighbors, n_jobs=10)
            knc.fit(X_log, labels)
            predictions = knc.predict(X_log)

            plt.figure(None, (16, 16))
            gs = plt.GridSpec(X.shape[1], X.shape[1])
            for i in range(X.shape[1]):
                for j in range(X.shape[1]):
                    if i == 1 and j == 1:
                        ax = plt.subplot(gs[i, j], adjustable='box', aspect=5.)
                        for c in range(np.max(labels) + 1):
                            plt.text(0.5, c / (np.max(labels) + 1), f"{c}".rjust(2),
                                     color="k", fontdict={"size": 5 + (40 / np.max(labels))},
                                     bbox=dict(facecolor=plt.cm.tab20(c / 19), edgecolor='none'),
                                     transform=ax.transAxes)
                        plt.axis("off")
                    elif i <= j:
                        pass
                    else:
                        ax = plt.subplot(gs[i, j])
                        sc = vcy.scatter_viz(X_log[ix, j], X_log[ix, i],
                                             s=7, alpha=0.2, lw=0, c=predictions[ix].astype(float), vmin=0, vmax=20, cmap=plt.cm.tab20)
                        plt.xlabel(feats_name[j])
                        plt.ylabel(feats_name[i])
                        for c in range(np.max(labels) + 1):
                            plt.text(np.median(X_log[ix, j][predictions[ix] == c]), np.median(X_log[ix, i][predictions[ix] == c]), f"{c}",
                                     color="k", fontdict={"size": 5})
            plt.tight_layout()
            plt.savefig(os.path.join(path_name, "QC_Class_Predicitions.png"), dpi=300)

            logging.info("Pickling the classifier object")
            pickle.dump(knc, open(os.path.join(path_name, "QC_Classifier.pickle"), "wb"))

        f = open(os.path.join(dm.paths().build, "QC_classifier", "names_qc_clusters.form.txt"), "w")
        f.write('\n'.join([f"{i:02}:" for i in range(np.max(labels) + 1)]))
        f.close()
