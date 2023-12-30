# This script is to create output files when the snakemake pipeline has not completely run.
# The output is created in the subdirectory ../temp/

#!/bin/bash

# creating list of vcfs
ls /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/*-biallelic.vcf.gz > ../temp/biallelic-vcfs.txt

# activating conda environment used in snakemake pipeline
source ~/.bashrc
conda activate /gpfs/project/projects/medbioinf/users/spani/scripts/1000g-ont/ont-1kg/snakemake_pipelines/augmented_graph/.snakemake/conda/91103c7d858528647eb601a6ea8dffb7_

# creating comma separate list of biallelic vcfs
csl_biallelic=`awk '{print $1}' ../temp/biallelic-vcfs.txt | paste -s -d, -`

# creating SV count plots
sv_count_script="/gpfs/project/projects/medbioinf/users/spani/scripts/1000g-ont/ont-1kg/snakemake_pipelines/augmented_graph/scripts/sv_count.py"
metadata="/gpfs/project/projects/medbioinf/users/spani/files/other/1000GP/igsr_sample_data.tsv"
python $sv_count_script -meta $metadata -vcf $csl_biallelic -output ../temp/ 2> ../temp/sv_count.log