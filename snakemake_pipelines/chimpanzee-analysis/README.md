# Description of the Pipeline

- Using minigraph, the chimp reference is aligned to the rGFA created in the graph augmentation part of the SAGA pipeline described.
- From the alignment, the alleles traversing the bubbles are extracted and then converted into a haploid VCF.
- The VCF is then processed into a BED file containing the information about the ancestral allele found in the chimp for each bubble.

## System Requirements

- python=3.12.2
- snakemake=7.32.4
- conda=24.7.1
- giggles (commit ID: [fd13a89682dbb3ebaabfa9eb68e6d3dc63638af4](https://github.com/samarendra-pani/giggles/tree/fd13a89682dbb3ebaabfa9eb68e6d3dc63638af4))
- gaftools (commit ID: [feaf7f45456fe49285fbfbd8e7cf1d2256dde65e](https://github.com/marschall-lab/gaftools/tree/feaf7f45456fe49285fbfbd8e7cf1d2256dde65e))

## Installation Guide

### Instructions

The pipeline is run with snakemake and the dependencies of each rule is resolved by the specified conda environment.

The following packages need installation:

- gaftools
    - link to github: https://github.com/marschall-lab/gaftools
    - for this pipeline, gaftools needs to be installed in a conda environment called `gaftools-env`. Instructions on conda installation is provided in the github page.

- giggles
    - link to github: https://github.com/samarendra-pani/giggles
    - for this pipeline, giggles needs to be installed in a conda environment called `giggles-env`. Instructions on conda installation is provided in the github wiki.
    - NOTE: This pipeline uses a later developmental commit of giggles which has the added subcommand of preparing phased panels. Users can download this development version and use it for the genotyping as well since this commit makes no changes to the genotyping algorithm.

## Instructions for use

### Input

Running the pipeline on our data will require downloading the entire dataset. Refer to the *config.yaml* file where it states the required data for the pipeline and where to download them. The required directories and files need to be specified in the variables and then the pipeline can be run.

### Output

The pipeline was run for various sets of synthetic reads produced from the chimpanzee reference assembly. The best result was obtained from a set of overlapping synthetic reads ranging from 4kbp to 500kbp.

The main outputs of this pipeline are:

- Annotated BED file with all the bubbles and its ancestral allele information at `results/4000-20000-50000-100000-200000-500000/ancestral-allele-annotations.bed`
- VCF with the alleles present in the chimp reference at `results/4000-20000-50000-100000-200000-500000/assembly-vcf/chimp.vcf`
- Alignment GAF for the chimp reference to the augmented rGFA at `results/4000-20000-50000-100000-200000-500000/assembly-to-graph-alignment/chimp.gaf`
- Updated callset VCFs. The phased and unphased callsets released as part of the study have been updated with the Ancestral Allele information
    - path to the updated phased VCF: `results/4000-20000-50000-100000-200000-500000/annotated-callset.phased.vcf`
    - path to the updated unphased VCF: `results/4000-20000-50000-100000-200000-500000/annotated-callset.unphased.vcf`

### Reproduction instructions

Under the directory `chimpanzee-analysis`, the command `snakemake --use-conda -j <number of cores> all` can be used to run the entire pipeline.

Due to the large requirements of certain jobs, a high performance computing architecture will be needed. For running the pipeline in such an architecture, refer to the snakemake documentation (https://snakemake.readthedocs.io/en/v7.32.4/) for further details.
