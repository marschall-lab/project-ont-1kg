# Description of the Pipeline

- Evaluating the SV recovery: how many SVs were recovered from the genotyping from the list of SVs introduced into the graph.
- The SAGA framework includes intermediate steps of introducing the discovered SVs into graph using pseudohaplotypes by minigraph. The consequent graph is then processed to identify the SV alleles. This pipeline disrupts the one-to-one identification of the SVs that were fed into the graph and the SVs that are retrieved after genotyping and phasing. Instead of trying to identify the SVs discovered by the various discovery softwares used, we identified the SV alleles that are exclusively found in the pseudohaplotypes constructed out of the discovered SVs. By extension, these pseudohaplotype-exclusive SVs act as a rough estimate for the discovered SVs. We compared the number of pseudohaplotype-exclusive SVs introduced into the genotyping panel and the fraction of them which were recovered after phasing.

## System Requirements

- python=3.12.2
- snakemake=7.32.4
- conda=24.7.1

## Installation Guide

### Instructions

The pipeline is run with snakemake and the dependencies of each rule is resolved by the specified conda environment.

## Instructions for use

### Input

Running the pipeline requires using files from the pipeline as defined in `config.yaml`. Some files are part of the release while others are intermediate files of the pipeline.

### Output

The main outputs of the pipeline are:

- Plot comparing SV counts of the pseudohaplotype-exclusive alleles present in the panel vs the callset at `results/counts.svg`
- Same comparison as an SV length distribution at `results/length-dist.svg`

### Reproduction instructions

Under the directory `sv-recovery`, the command `snakemake --use-conda all` can be used to run the entire pipeline.