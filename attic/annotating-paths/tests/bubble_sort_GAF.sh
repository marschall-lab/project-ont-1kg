#!/bin/bash
#PBS -A LongReadPanGenie
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -l walltime=01:00:00

# source bash
source $HOME/.bashrc
source $HOME/.bash_profile
conda activate giggles

python /gpfs/project/projects/medbioinf/users/spani/scripts/1000g-ont/ont-1kg/snakemake_pipelines/annotating-paths/scripts/bubble-sort.py --gfa /gpfs/project/projects/medbioinf/users/spani/results/1000GP/annotating-paths/bubble_calling_and_tagging/chm13-90c.r518_tagged_noseq.gfa --gaf /gpfs/project/projects/medbioinf/users/spani/NA12878.gaf --output /gpfs/project/projects/medbioinf/users/spani/NA12878_sorted.gaf 2> /gpfs/project/projects/medbioinf/users/spani/NA12878_sorted.log
