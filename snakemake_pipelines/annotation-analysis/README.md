# Description of the Pipeline

- Get some statistics from the genotypes after annotated by SVAN

## System Requirements

- python=3.12.2
- snakemake=7.32.4
- conda=24.7.1

## Installation Guide

### Instructions

The pipeline is run with snakemake and the dependencies of each rule is resolved by the specified conda environment.

## Instructions for use

### Input

Refer to the *config.yaml* file where it states the required data for the pipeline and where to download them. The required directories and files need to be specified in the variables and then the pipeline can be run.

### Output

It produces the total count of NUMT, PSD, and INV_DUP sites per sample.

### Reproduction instructions

Under the directory `annotation-analysis`, the command `snakemake --use-conda -j <number of cores> all` can be used to run the entire pipeline.

Due to the large requirements of certain jobs, a high performance computing architecture will be needed. For running the pipeline in such an architecture, refer to the snakemake documentation (https://snakemake.readthedocs.io/en/v7.32.4/) for further details.