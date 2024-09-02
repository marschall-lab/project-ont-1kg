# Description of the Pipeline

- Doing the Vamos Analysis for VNTRs using our aligned BAMs with the sites list from Gustafson, Gibson, Damaraju et al, 2024.
- Generated Vamos results for the ONT data and compared the two VCFs with the overlapping samples (HG03499, HG03548).

## System Requirements

- python=3.12.2
- snakemake=7.32.4
- conda=24.7.1

## Installation Guide

### Instructions

The pipeline is run with snakemake and the dependencies of each rule is resolved by the specified conda environment.

## Instructions for use

### Input

Running the pipeline on our data will require downloading the CRAM files for the two samples in common with Gustafson, Gibson, Damaraju et al, 2024 (the common samples are mentioned in the config file). Refer to the *config.yaml* file where it states the required data for the pipeline and where to download them. The required directories and files need to be specified in the variables and then the pipeline can be run.

### Output

The main outputs of the pipeline are:

- Plots using an allele distance metric as a histogram comparing the two vcfs with the two samples.

### Reproduction instructions

Under the directory `chimpanzee-analysis`, the command `snakemake --use-conda -j <number of cores> all` can be used to run the entire pipeline.

Due to the large requirements of certain jobs, a high performance computing architecture will be needed. For running the pipeline in such an architecture, refer to the snakemake documentation (https://snakemake.readthedocs.io/en/v7.32.4/) for further details.