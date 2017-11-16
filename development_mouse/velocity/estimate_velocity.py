from typing import *
import os
import loompy
from scipy import sparse
import numpy as np
import cytograph as cg
import development_mouse as dm
import velocyto as vcy
import matplotlib.pyplot as plt
import luigi


class EstimateVelocity(luigi.Task):
    """Luigi Task to run velocyto
    """
    tissue = luigi.Parameter(description="name of the tissue from tool specification file.",
                             always_in_help=True)

    def requires(self) -> List[luigi.Task]:
        """
        Arguments
        ---------
        `AggregateL1`:
            passing ``tissue``
        `ClusterL1`:
            passing ``tissue``
        """
        # NOTE: not sure it needs AggregateL1
        return [dm.ClusterL1(tissue=self.tissue),
                dm.AggregateL1(tissue=self.tissue)]

    def output(self) -> luigi.Target:
        """
        Returns
        -------
        folder: ``L1_[TISSUE]_velocity``:
            Note this is kind of a hack to luigi, single files will not be regenerated but whole folder will.
        """
        return luigi.LocalTarget(os.path.join(dm.paths().build, "L1_" + self.tissue + "_velocity"))

    def run(self) -> None:
        """
        Reads the output of `AggregateL1` and runs velocyto analysis:
        """
        logging = cg.logging(self, True)
        logging.info("Exporting cluster data")
        with self.output().temporary_path() as out_dir:
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            
            logging.info("Loading loom file in memory as a VelocytoLoom object")
            vlm = vcy.VelocytoLoom(self.input()[0].fn)
            vlm.set_clusters(cluster_labels=vlm.ca["ClusterName"])
            logging.info("Plotting report on spliced, ambiguous, unpliced fraction")
            vlm.plot_fractions(save2file=os.path.join(out_dir, "L1_" + self.tissue + "_sau_fractions.pdf"))

            # NOTE: code below is basically identical to `default_filter_and_norm` but with the exception of adjust_totS_totU
            
            # Heuristics, we should set better heuristic and could add a config file with parameters for analysis
            max_expr_avg = 40
            min_expr_counts = max(20, min(100, vlm.S.shape[1] * 2.25e-3))
            min_cells_express = max(10, min(50, vlm.S.shape[1] * 1.5e-3))
            N = max(1000, min(int((vlm.S.shape[1] / 1000)**(1 / 3) / 0.0008), 5000))
            min_avg_U = 0.01
            min_avg_S = 0.08

            # NOTE: not sure if this is needed with the new init
            vlm.normalize("S", size=True, log=False)
            vlm.normalize("U", size=True, log=False)

            logging.info("Performing gene filtering by S detection")
            vlm.score_detection_levels(min_expr_counts=min_expr_counts, min_cells_express=min_cells_express)
            vlm.filter_genes(by_detection_levels=True)

            logging.info("Performing gene filtering by Cv vs mean relation")
            vlm.score_cv_vs_mean(N=N, max_expr_avg=max_expr_avg)
            vlm.filter_genes(by_cv_vs_mean=True)

            logging.info("Performing gene filtering by U detection")
            vlm.score_detection_levels(min_expr_counts=0, min_cells_express=0,
                                       min_expr_counts_U=int(min_expr_counts / 2) + 1,
                                       min_cells_express_U=int(min_cells_express / 2) + 1)
            
            if hasattr(vlm, "cluster_labels"):
                logging.info("Performing gene filtering by cluster expression")
                vlm.score_cluster_expression(min_avg_U=min_avg_U, min_avg_S=min_avg_S)
                vlm.filter_genes(by_detection_levels=True, by_cluster_expression=True)
            else:
                vlm.filter_genes(by_detection_levels=True)

            vlm.normalize_by_total(plot=True)

            logging.info("Preparing dataset for velocity extimation")
            vlm.perform_PCA()
            n_comps = int(np.where(np.diff(np.diff(np.cumsum(vlm.pca.explained_variance_ratio_)) > 0.002))[0][0])
            n_comps = min(n_comps, 50)
            k = int(min(1000, max(10, np.ceil(vlm.S.shape[1] * 0.02))))

            logging.info(f"Considering {n_comps} components and {k} nearest neighbours")
            vlm.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k * 8, b_maxl=k * 4, n_jobs=8)

            vlm.normalize_median()  # NOTE: it had problems in the past...

            logging.info(f"Fitting gammas for {vlm.Sx_sz.shape[1]} genes")
            vlm.fit_gammas()

            logging.info("Calculate velocity")
            vlm.predict_U()
            vlm.calculate_velocity()
            vlm.calculate_shift(assumption="constant_velocity")
            vlm.extrapolate_cell_at_t(delta_t=1)  # NOTE: we should determine delta t in a better way

            vlm.ts = np.column_stack([vlm.ca["_X"], vlm.ca["_Y"]])  # load the embedding from previous analysis

            logging.info(f"Estimating transition probability, this step will require time")
            n_neighbors = int(vlm.S.shape[1] / 5)
            vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt",
                                         n_neighbors=n_neighbors, knn_random=True, sampled_fraction=1, n_jobs=16)
            # NOTE here we might want to tune the number of jobs used in relation to the number of concurrent tasks
            # NOTE here we might want to change the sampled fraction to a lower number to make things faster

            vlm.calculate_embedding_shift(sigma_corr=0.05)  # NOTE: this parameter could be tuned

            vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=300)  # NOTE: this parameters could be tuned

            plt.figure(None, (20, 20))
            vlm.plot_grid_arrows(scatter_kwargs_dict={"alpha": 0.35, "lw": 0.35, "edgecolor": "0.4", "s": 38, "rasterized": True},
                                 min_mass=10, angles='xy', scale_units='xy',
                                 headaxislength=2.75, headlength=5, headwidth=4.8, quiver_scale=0.25)
            # NOTE: this parameters could be tuned. In particular min_mass!
            plt.savefig(os.path.join(out_dir, "L1_" + self.tissue + "_velocity_TSNE.png"))
