#!/bin/bash
#PBS -A LongReadPanGenie
#PBS -l select=1:ncpus=1:mem=170gb
#PBS -l walltime=71:00:00

# creating list of vcfs
ls /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/*-multiallelic.vcf.gz > ../temp/multiallelic-vcfs.txt

# creating merged multiallelic vcf
source ~/.bashrc
conda activate base
bcftools merge --no-version -l ../temp/multiallelic-vcfs.txt -o /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/multisample-multiallelic.vcf

# gziping
conda activate /gpfs/project/projects/medbioinf/users/spani/scripts/1000g-ont/ont-1kg/snakemake_pipelines/augmented_graph/.snakemake/conda/91103c7d858528647eb601a6ea8dffb7_
bgzip -c /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/multisample-multiallelic.vcf > /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/multisample-multiallelic.vcf.gz
tabix -p vcf /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/multisample-multiallelic.vcf.gz

# removing list of files
rm ../temp/multiallelic-vcfs.txt

# converting to biallelic
cat /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/multisample-multiallelic.vcf | python convert-to-biallelic.py /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/panel/giggles-ready_biallelic.vcf.gz | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n "}' > /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/multisample-biallelic.vcf
bgzip -c /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/multisample-biallelic.vcf > /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/multisample-biallelic.vcf.gz
tabix -p vcf /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/multisample-biallelic.vcf.gz

# removing ungzipped files
rm /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/multisample-multiallelic.vcf
rm /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/multisample-biallelic.vcf
