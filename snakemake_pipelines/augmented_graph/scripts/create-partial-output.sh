# This script is to create output files when the snakemake pipeline has not completely run.
# The output is created in the subdirectory ../temp/

#!/bin/bash

# creating list of vcfs
ls /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/*-biallelic.vcf.gz > ../temp/biallelic-vcfs.txt
ls /gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/genotypes/*-multiallelic.vcf.gz > ../temp/multiallelic-vcfs.txt

# activating conda environment used in snakemake pipeline
conda activate /gpfs/project/projects/medbioinf/users/spani/scripts/1000g-ont/ont-1kg/snakemake_pipelines/augmented_graph/.snakemake/conda/91103c7d858528647eb601a6ea8dffb7_

# creating merged multiallelic vcf
bcftools merge --no-version -l ../temp/multiallelic-vcfs.txt -o ../temp/multisample-multiallelic.vcf
bgzip -c ../temp/multisample-multiallelic.vcf > ../temp/multisample-multiallelic.vcf.gz
tabix -p vcf ../temp/multisample-multiallelic.vcf.gz

# get the biallelic multisample VCF
convert_to_biallelic="/gpfs/project/projects/medbioinf/users/spani/scripts/1000g-ont/ont-1kg/snakemake_pipelines/augmented_graph/scripts/convert-to-biallelic.py"
biallelic_panel="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/panel/giggles-ready_biallelic.vcf.gz"
cat ../temp/multisample-multiallelic.vcf | python $convert_to_biallelic $biallelic_panel | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k2,2n \"}' > ../temp/multisample-biallelic.vcf
bgzip -c ../temp/multisample-biallelic.vcf > ../temp/multisample-biallelic.vcf.gz
tabix -p vcf ../temp/multisample-biallelic.vcf.gz

# creating comma separate list of biallelic vcfs
csl_biallelic=`awk '{print $1}' ../temp/biallelic-vcfs.txt | paste -s -d, -`

# creating SV count plots
sv_count_script="/gpfs/project/projects/medbioinf/users/spani/scripts/1000g-ont/ont-1kg/snakemake_pipelines/augmented_graph/scripts/sv_count.py"
metadata="/gpfs/project/projects/medbioinf/users/spani/files/other/1000GP/igsr_sample_data.tsv"
python $sv_count_script -meta $metadata -vcf $csl_biallelic -output ../temp/ 2> ../temp/sv_count.log

# collecting vcf stats
script="/gpfs/project/projects/medbioinf/users/spani/scripts/1000g-ont/ont-1kg/snakemake_pipelines/augmented_graph/scripts/collect-vcf-stats.py"
metadata="/gpfs/project/projects/medbioinf/users/spani/files/other/1000GP/igsr_sample_data.tsv"
biallelic_panel="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/panel/giggles-ready_biallelic.vcf.gz"
callset="../temp/multisample-biallelic.vcf.gz"
python $script -meta $metadata -panel $biallelic_panel -callset $callset > ../temp/callset-stats.tsv 2> ../temp/callset-stats.log

# adding bubble ids to stats
script="/gpfs/project/projects/medbioinf/users/spani/scripts/1000g-ont/ont-1kg/snakemake_pipelines/augmented_graph/scripts/add-bub-info.py"
table="../temp/callset-stats.tsv"
multiallelic_panel="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/minigraph-extended_all/panel/giggles-ready_multiallelic.vcf.gz"
python $script -table $table -panel $multiallelic_panel -output ../temp/variant-stats.tsv

# plotting HWE and AF
script="/gpfs/project/projects/medbioinf/users/spani/scripts/1000g-ont/ont-1kg/snakemake_pipelines/augmented_graph/scripts/plot-vcf-stats.py"
python $script -table ../temp/variant-stats.tsv -output ../temp/