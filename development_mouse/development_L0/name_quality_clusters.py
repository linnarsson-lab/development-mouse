import luigi
import os
import development_mouse as dm


class NameQualityClusters(luigi.ExternalTask):
    """External task to allow manual naming of clusters found by MakeQualityClassifier

    The user is meant to run MakeQualityClassifier first, this will generate ``names_qc_clusters.form.txt``.
    This file needs to be modified to contain cluster-number:clustername pairs and renamed ``names_qc_clusters.txt``
    An example of the ``names_qc_clusters.txt`` is illusterated below:

    Example
    -------
    ```
    0:HQ1
    1:HQ2
    2:LQ1
    3:LQ2
    4:DE1
    5:NU1
    6:RBC
    ```

    """
    def output(self) -> luigi.Target:
        return luigi.LocalTarget(os.path.join(dm.paths().build, "names_qc_clusters.txt"))
