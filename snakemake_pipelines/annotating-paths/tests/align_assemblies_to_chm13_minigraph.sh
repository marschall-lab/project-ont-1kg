#!/bin/bash
#PBS -A LongReadPanGenie
#PBS -l select=1:ncpus=24:mem=100gb
#PBS -l walltime=10:59:00
#PBS -o /gpfs/project/projects/medbioinf/users/spani/results/1000GP/annotating-paths/tests/$PBS_JOBNAME.out
#PBS -e /gpfs/project/projects/medbioinf/users/spani/results/1000GP/annotating-paths/tests/$PBS_JOBNAME.err 

/gpfs/project/projects/medbioinf/users/spani/packages/minigraph/minigraph --vc -cx lr /gpfs/project/projects/medbioinf/users/spani/files/gfa/1000GP/chm13-90c.r518.gfa.gz /gpfs/project/projects/medbioinf/data/hprc/assemblies/HG01258.maternal.f1_assembly_v1.fa.gz -t 16 > /gpfs/project/projects/medbioinf/users/spani/results/1000GP/annotating-paths/tests/HG01258.maternal.gaf
