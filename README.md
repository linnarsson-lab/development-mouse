# development-mouse
Analysis pipeline for the mouse developing nervous system project.

Publication: [Braun et al: Comprehensive cell atlas of the first-trimester developing human brain](https://www.science.org/doi/10.1126/science.adf1226)

Download data:
  * Complete dataset, 31053 genes x 292495 cells, [anndata](https://anndata.scverse.org/en/stable/) format: [developing_mouse_nervous_system.h5ad](https://www.dropbox.com/scl/fi/taxtkhr3fezv6qoe86vcb/developing_mouse_nervous_system.h5ad?rlkey=fnphj9kokc2mstzudx3cmze6n&st=u7sawgjl&dl=0)
  * One column per cluster, 31053 genes x 798 clusters, [loom](https://loompy.org) format: [dev_all.agg.loom](https://www.dropbox.com/scl/fi/67mu8p6sxqymui39ltz54/dev_all.agg.loom?rlkey=64invw7it2hqe578k1jidpztt&st=rp3jocip&dl=0)
  * Spatial transcriptomics at 5w: [640um](https://www.dropbox.com/scl/fi/4igq3gd8g6qj5nqo2rvdh/LBEXP20211117_EEL_HE_5w_640um_RNA_transformed_assigned.parquet?rlkey=y0vrnb9c05w9wjjpcq52ophjo&st=sblutml4&dl=0)
 [810um](https://www.dropbox.com/scl/fi/obtfabltd9b4fd1r27tbi/LBEXP20211119_EEL_HE_5w_810um_RNA_transformed_assigned.parquet?rlkey=ot9y67j85wpw0988eueh35uy0&st=r5h9r66a&dl=0)
 [970um](https://www.dropbox.com/scl/fi/66kucqzt1x4mwp65q1l85/LBEXP20211113_EEL_HE_5w_970um_RNA_transformed_assigned.parquet?rlkey=ltlss5qz1ef188dvahnn8qvf0&st=a49q7nlp&dl=0)

## Installation

1. Install [luigi](https://luigi.readthedocs.io/en/stable/), [cytograph](https://github.com/linnarsson-lab/cytograph) and [loompy](https://github.com/linnarsson-lab/loompy)

2. Clone the repository to your computer:

```
git clone https://github.com/linnarsson-lab/development-mouse.git
```

3. Install in development mode using `pip`:

```
cd development-mouse
pip install -e .
```

## Setting up your environment

In order to run the pipeline, you will need to set up your samples folder:

* A folder containing the raw samples (loom files named like `10X23_1.loom`)
* Inside it, a sub-folder named `classified` and containing pre-classified clusters (loom files named like `L0_Cortex1.loom`)
* Also inside it, a sub-folder named `metadata` and containing a file named `metadata.xlsx` with sample metadata

For example, you could have something like this:

```
loom_samples/
    10X23_1.loom
    10X23_2.loom
    10X53_1.loom
    ...
    classified/
        L0_Cortex1.loom
        L0_Hippocampus.loom
        ...
        classifier.pickle
    metadata/
        metadata.xlsx
```

**Note:** On monod, a samples folder is available at `/data/proj/chromium/`. Please use this directly instead of making a copy.

**Note:** Instead of sample files directly in the top-level folder, sample files can also be in subfolders named after the sample ID, or in a sub-folder of such a subfolder, named `velocyto`. Thus, these are all equivalent:

```
10X43_2.loom
10X43_2/10X43_2.loom
10X43_2/velocyto/10X43_2.loom
```

If it exists, the sample in `velocyto` will take precedence.

Furthermore, you need a file in the current directory (typically, `development-mouse`) named `pooling_specification.tab`, which gives a list of all the samples and their pool names. This file has four columns: *SampleID* (like `10X04_1`), *Pool* (like `Hippocampus`), *TimepointPool* (always `none`), *QC* (`OK` or `FAILED`), and *Project* (`Adolescent`). 

Samples with `QC == FAILED` will be ignored for all analyses.


## Running the pipeline

1. Create a folder to hold the output of the build, e.g. `mkdir /data/proj/development/build_20171107`.

2. Run `luigi`. For example:

```
luigi --workers 15 --local-scheduler --module development_mouse TASKNAME --paths-samples /data/proj/chromium/loom_samples/ --paths-build /data/proj/development/build_20171107
```

The command line arguments are as follows:

Argument|Effect
----|----
`--workers` | Number of parallel process to use
`--local-scheduler` |run the pipeline directly, not in client/server mode.
`--module development_mouse TASKNAME`| run the `TASKNAME` task in the Python module `development_mouse`.
`--paths-samples /data/proj/chromium/loom_samples/` | a configuration of the `paths` object, setting the sample path
`--paths-build /data/proj/development/build_20171107` | a configuration of the `paths` object, setting the build path

Tasks Descriptions
==================

# Level1 tasks

[TODO] Add image of the tasks hierarchy

## Level1Analysis

It is the main task to run. It triggers ExportL1 that in turn will trigger all level 1 set.

# Velocity related

below an example how to run velocity on a single sample

## Preliminary Level1 Analysis

If we want to make make sure that Level1 Analysis has been appropriatelly run, run:

```
luigi --local-scheduler --module development_mouse ExportL1 --tissue Forebrain_E9-11 --paths-samples /data/proj/chromium/loom_samples/ --paths-build /data/proj/development/build_20171115
```

## EstimateVelocity

```
nohup luigi --local-scheduler --module development_mouse EstimateVelocity --tissue Forebrain_E9-11 --paths-samples /data/proj/chromium/ --paths-build /data/proj/development/build_20171115 > ~/luigi_EstimateVelocity_EstimateVelocity.log &
```
