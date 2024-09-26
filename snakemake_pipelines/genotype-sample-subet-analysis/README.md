# Description of the Pipeline

- Breaking down the sample list into different coverage ranges, and n50 ranges to investigate the effect of coverage on results
- The pipeline looks at the following
    - The quality of the genotypes (before the population phasing step by ShapeIt5)
    - Difference in SVs discovered vs genotyped
    - [TODO] Downsampling and genotyping of high coverage samples to show effects of downsampling.

## System Requirements

- python=3.12.2
- snakemake=7.32.4
- conda=24.7.1

## Installation Guide

### Instructions

The pipeline is run with snakemake and the dependencies of each rule is resolved by the specified conda environment.

## Instructions for use

### Input

Running the pipeline on our data will creating and downloading the data specified in the `config.yaml` file. Some data are available as part of the data release for this study while others are intermediate files of the pipeline and data downloadable from public sources.

### Output

The main outputs of this pipeline are:

- QC plots for various coverage ranges at `results/coverage-experiments/qc-plots/all.pdf`
- SV counts for discovery vs genotyping plot at `results/coverage-experiments/sv_count_per_sample/plot.svg`
- QC plots for various n50 ranges at `results/n50-experiments/qc-plots/all.pdf`
- SV counts for discovery vs genotyping plot at `results/n50-experiments/sv_count_per_sample/plot.svg`

### Reproduction instructions

Under the directory `genotype-sample-subset-analysis`, the command `snakemake --use-conda` can be used to run the entire pipeline.
