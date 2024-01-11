# This script is to create output files when the snakemake pipeline has not completely run.
# The output is created in the subdirectory ../temp/

#!/bin/bash

# activating conda environment used in snakemake pipeline
source ~/.bashrc
conda activate /gpfs/project/projects/medbioinf/users/spani/scripts/1000g-ont/ont-1kg/snakemake_pipelines/augmented_graph/.snakemake/conda/91103c7d858528647eb601a6ea8dffb7_

# creating comma separate list of biallelic vcfs
csl_biallelic=`awk '{print $1}' biallelic-vcfs.txt | paste -s -d, -`

# creating SV count plots
sv_count_script="/gpfs/project/projects/medbioinf/users/spani/scripts/1000g-ont/ont-1kg/snakemake_pipelines/svarp-giggles_combined_pipeline/scripts/sv_length.py"
metadata="/gpfs/project/projects/medbioinf/users/spani/files/other/1000GP/igsr_sample_data.tsv"
python $sv_count_script -meta $metadata -vcf $csl_biallelic -output /gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/