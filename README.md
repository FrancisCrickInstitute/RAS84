
# RAS84 Analysis Code

This Github repo contains the code needed to download and preprocess
the data, build the SVM classifier and regenerate all the figures
associated with the publication [RAS oncogenic activity predicts
response to chemotherapy and outcome in lung
adenocarcinoma](https://www.nature.com/articles/s41467-022-33290-0)

## Running the code

First step, clone the repo.

```

git clone https://github.com/FrancisCrickInstitute/RAS84.git
cd RAS84

```

To run the analysis code you will need `R-3.6.0` installed along with
the libraries contained in `renv.lock`. The `renv` environment can be
[reinitialised](https://rstudio.github.io/renv/articles/collaborating.html)
by running `R` within the `RAS84` directory and running `renv::restore()`


## The Data

`data/scripts` contains shell scripts to download each of the
dependent datasets. The paths are configured to run these scripts from
within `data/scripts`. Raw data files are downloaded to
`data/downloads`

There are also R scripts to preprocess the raw data and build
`SummarizedExperiment` objects. The SE objects are written to
`data/objects`

## The Analysis Scripts

The following analysis `Rmd` scripts can be found in the top level of
the repo.

1. 1_build_RAS84_CCLE.Rmd - The CCLE analysis to build RAS84 from the
parent signatures

2. 2_drug_screen_GDSC_CTRP_CCLE.Rmd - The CCLE drug screen analysis

3. 3_RAS84_patient_classification_TCGA_LUAD.Rmd - TCGA LUAD classification analysis

4. 4_RAS84_RAG_mutation_analysis_TCGA_LUAD.Rmd - TCGA LUAD RAS Activity Group (RAG) mutation analysis

5. 5_survival_analysis_TCGA_LUAD.Rmd - TCGA LUAD and Uppsala survival
   analysis
   
6. 6_RAS84_pancancer_analysis_TCGA.Rmd - TCGA pancancer RAS84 analysis

7. 7_RAS84_pancancer_survival_analysis.Rmd - TCGA pancancer RAS84
survival analysis

8. build_SVM.R - SVM build script

All the `Rmd` scripts contain global variables in the first code block
that set the paths to the required data resources and output paths for the
figures. All the defined data resources must be downloaded and build
using the download and init scripts described above before the
analysis scripts can be run. All figures are written to `figures/`. 


