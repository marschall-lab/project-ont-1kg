.ONESHELL:

# Need to specify bash in order for conda activate to work.
SHELL=/bin/bash
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

pre-filtering: plot-hwe-af-unfiltered boxplots-unfiltered rausch-curves-unfiltered audano-curves-unfiltered

tag-filter-vcfs:
	bcftools view --min-ac 1 original-graph/multisample-genotypes/multisample-biallelic.vcf.gz | bcftools +fill-tags | bgzip -c > original-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz
	bcftools view --min-ac 1 augmented-graph/multisample-genotypes/multisample-biallelic.vcf.gz | bcftools +fill-tags | bgzip -c > augmented-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz
	bcftools view -S unrel.samples.tsv original-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz | bgzip -c > original-graph/multisample-genotypes/biallelic.ac0-filtered.unrelated-samples.vcf.gz
	bcftools view -S unrel.samples.tsv augmented-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz | bgzip -c > augmented-graph/multisample-genotypes/biallelic.ac0-filtered.unrelated-samples.vcf.gz

map-bub-ids:
	python scripts/map-bub-ids-to-allele-ids.py -vcf original-graph/multisample-genotypes/multisample-multiallelic.vcf.gz > original-graph/map-bub_ids-to-allele_ids.tsv
	python scripts/map-bub-ids-to-allele-ids.py -vcf augmented-graph/multisample-genotypes/multisample-multiallelic.vcf.gz > augmented-graph/map-bub_ids-to-allele_ids.tsv

boxplots-unfiltered:
	python scripts/boxplots-lengths.py -meta igsr_sample_data.tsv -coverages 1k_ont_data_overview-samples_and_sequencing.tsv -names original,augment -vcfs original-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz,augmented-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz -output joint-plots/unfiltered/boxplot-lengths/
	python scripts/boxplots-counts.py -meta igsr_sample_data.tsv -coverages 1k_ont_data_overview-samples_and_sequencing.tsv -names original,augment -vcfs original-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz,augmented-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz -output joint-plots/unfiltered/boxplot-counts/

rausch-curves-unfiltered:
	python scripts/rausch-curves.py -sample-sheet sample_sheet.tsv -names original,augment -vcfs original-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz,augmented-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz -output joint-plots/unfiltered/rausch-curves/

audano-curves-unfiltered:
	python scripts/audano-growth-curves.py -vcf original-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz -meta igsr_sample_data.tsv -map original-graph/map-bub_ids-to-allele_ids.tsv -sample-sheet sample_sheet.tsv -output original-graph/plots/growth-curves
	python scripts/audano-growth-curves.py -vcf augmented-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz -meta igsr_sample_data.tsv -map augmented-graph/map-bub_ids-to-allele_ids.tsv -sample-sheet sample_sheet.tsv -output augmented-graph/plots/growth-curves

create-sample-list-table:
	python scripts/sample-list-supplementary-table.py > supplementary_table-sample_list.tsv


# create-hwe-af-popwise-stats:
# 	python scripts/hwe-af-popwise-stats.py -sample-sheet sample_sheet.tsv -meta igsr_sample_data.tsv -bi-panel augmented-graph/prepared-panel/giggles-ready_biallelic.vcf.gz -multi-panel augmented-graph/prepared-panel/giggles-ready_multiallelic.vcf.gz -callset augmented-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz -coverages 1k_ont_data_overview-samples_and_sequencing.tsv -output hwe-af-stats/augmented
# 	python scripts/hwe-af-popwise-stats.py -sample-sheet sample_sheet.tsv -meta igsr_sample_data.tsv -bi-panel original-graph/prepared-panel/panel-biallelic.vcf.gz -multi-panel original-graph/prepared-panel/panel-multiallelic.vcf.gz -callset original-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz -coverages 1k_ont_data_overview-samples_and_sequencing.tsv -output hwe-af-stats/original

# plot-hwe-af-popwise-stats:
# 	python scripts/hwe-af-popwise-plot.py -table hwe-af-stats/augmented-all_cov.tsv -output plots/augmented/hwe_af-all
# 	python scripts/hwe-af-popwise-plot.py -table hwe-af-stats/augmented-low_cov.tsv -output plots/augmented/hwe_af-low
# 	python scripts/hwe-af-popwise-plot.py -table hwe-af-stats/augmented-high_cov.tsv -output plots/augmented/hwe_af-high
# 	python scripts/hwe-af-popwise-plot.py -table hwe-af-stats/original-all_cov.tsv -output plots/original/hwe_af-all
# 	python scripts/hwe-af-popwise-plot.py -table hwe-af-stats/original-low_cov.tsv -output plots/original/hwe_af-low
# 	python scripts/hwe-af-popwise-plot.py -table hwe-af-stats/original-high_cov.tsv -output plots/original/hwe_af-high

qc-table-unfiltered:
	python scripts/qc-table-creation.py -sample-sheet sample_sheet.tsv -bi-panel augmented-graph/prepared-panel/giggles-ready_biallelic.vcf.gz -multi-panel augmented-graph/prepared-panel/giggles-ready_multiallelic.vcf.gz -bi-callset augmented-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz > augmented-graph/qc-table.tsv
	python scripts/qc-table-creation.py -sample-sheet sample_sheet.tsv -bi-panel original-graph/prepared-panel/panel-biallelic.vcf.gz -multi-panel original-graph/prepared-panel/panel-multiallelic.vcf.gz -bi-callset original-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz > original-graph/qc.final.tsv
	
merge-tables-unfiltered:
	python scripts/merge-qc-tables.py augmented-graph/qc-table.tsv augmented-graph/qc-phased.table.tsv augmented-graph/qc.final.tsv
	rm augmented-graph/qc-table.tsv

plot-hwe-af-unfiltered:
	python scripts/hwe-af-simple.py -table original-graph/qc.final.tsv -output original-graph/plots/hwe-af/
	python scripts/hwe-af-simple.py -table augmented-graph/qc.final.tsv -output augmented-graph/plots/hwe-af/


### BCFtools based filtering

bcftools-based-regression: bcftools-filtering bcftools-filtered-qc-table bcftools-merge-table

# cannot run directly here since it does not support conda activate
bcftools-filtering:
	$(CONDA_ACTIVATE) bcftools-latest
	bcftools +fill-tags augmented-graph/multisample-genotypes/multisample-biallelic.vcf.gz -- -t all | bcftools filter -O z -o genotype-filtering/augmented/biallelic.maskgq10.dropsites.vcf.gz -i 'INFO/AC >= 1 && INFO/F_MISSING <= 0.05' -
	bcftools view --min-ac 1 -m2 -M2 genotype-filtering/augmented/biallelic.maskgq10.dropsites.vcf.gz  | bcftools +setGT - -- -t q -n . -i 'FMT/GQ<10' | bcftools +fill-tags - -- -t all | bcftools filter -O z -o genotype-filtering/augmented/biallelic.filtered.all-samples.vcf.gz -i 'INFO/AC >= 1 && INFO/F_MISSING <= 0.05' -
	rm genotype-filtering/augmented/biallelic.maskgq10.dropsites.vcf.gz
	bcftools view -S unrel.samples.tsv genotype-filtering/augmented/biallelic.filtered.all-samples.vcf.gz | bgzip -c > genotype-filtering/augmented/biallelic.filtered.unrel-samples.vcf.gz
	bcftools +fill-tags original-graph/multisample-genotypes/multisample-biallelic.vcf.gz -- -t all | bcftools filter -O z -o genotype-filtering/original/biallelic.maskgq10.dropsites.vcf.gz -i 'INFO/AC >= 1 && INFO/F_MISSING <= 0.05'
	bcftools view --min-ac 1 -m2 -M2 genotype-filtering/original/biallelic.maskgq10.dropsites.vcf.gz  | bcftools +setGT - -- -t q -n . -i 'FMT/GQ<10' | bcftools +fill-tags - -- -t all | bcftools filter -O z -o genotype-filtering/original/biallelic.filtered.all-samples.vcf.gz -i 'INFO/AC >= 1 && INFO/F_MISSING <= 0.05' -
	rm genotype-filtering/original/biallelic.maskgq10.dropsites.vcf.gz
	bcftools view -S unrel.samples.tsv genotype-filtering/original/biallelic.filtered.all-samples.vcf.gz | bgzip -c > genotype-filtering/original/biallelic.filtered.unrel-samples.vcf.gz
	$(CONDA_ACTIVATE) base

# we need to generate the qc table for the other codes
bcftools-filtered-qc-table:
	python scripts/qc-table-creation.py -sample-sheet sample_sheet.tsv -bi-panel augmented-graph/prepared-panel/giggles-ready_biallelic.vcf.gz -multi-panel augmented-graph/prepared-panel/giggles-ready_multiallelic.vcf.gz -bi-callset genotype-filtering/augmented/biallelic.filtered.all-samples.vcf.gz > genotype-filtering/augmented/qc-table.tsv
	python scripts/qc-table-creation.py -sample-sheet sample_sheet.tsv -bi-panel original-graph/prepared-panel/panel-biallelic.vcf.gz -multi-panel original-graph/prepared-panel/panel-multiallelic.vcf.gz -bi-callset genotype-filtering/original/biallelic.filtered.all-samples.vcf.gz > genotype-filtering/original/qc.final.tsv

# need to merge the table with the stats from 
bcftools-merge-table:
	python scripts/merge-qc-tables.py genotype-filtering/augmented/qc-table.tsv augmented-graph/qc-phased.table.tsv genotype-filtering/augmented/qc.final.tsv
	rm genotype-filtering/augmented/qc-table.tsv



### All the plotting scripts after the filtering is done

post-filtering: plot-hwe-af-filtered boxplots-filtered rausch-curves-filtered audano-curves-filtered boxplots-all rausch-curves-all

plot-hwe-af-filtered:
	python scripts/hwe-af-simple.py -table genotype-filtering/augmented/qc.final.tsv -output genotype-filtering/augmented/plots/hwe-af/
	python scripts/hwe-af-simple.py -table genotype-filtering/original/qc.final.tsv -output genotype-filtering/original/plots/hwe-af/

boxplots-filtered:
	python scripts/boxplots-lengths.py -meta igsr_sample_data.tsv -coverages 1k_ont_data_overview-samples_and_sequencing.tsv -names original-filtered,augment-filtered -vcfs genotype-filtering/original/biallelic.filtered.all-samples.vcf.gz,genotype-filtering/augmented/biallelic.filtered.all-samples.vcf.gz -output joint-plots/filtered/boxplot-lengths/
	python scripts/boxplots-counts.py -meta igsr_sample_data.tsv -coverages 1k_ont_data_overview-samples_and_sequencing.tsv -names original-filtered,augment-filtered -vcfs genotype-filtering/original/biallelic.filtered.all-samples.vcf.gz,genotype-filtering/augmented/biallelic.filtered.all-samples.vcf.gz -output joint-plots/filtered/boxplot-counts/

rausch-curves-filtered:
	python scripts/rausch-curves.py -sample-sheet sample_sheet.tsv -names original-filtered,augment-filtered -vcfs genotype-filtering/original/biallelic.filtered.all-samples.vcf.gz,genotype-filtering/augmented/biallelic.filtered.all-samples.vcf.gz -output joint-plots/filtered/rausch-curves/

audano-curves-filtered:
	python scripts/audano-growth-curves.py -vcf genotype-filtering/original/biallelic.filtered.all-samples.vcf.gz -meta igsr_sample_data.tsv -map original-graph/map-bub_ids-to-allele_ids.tsv -sample-sheet sample_sheet.tsv -output genotype-filtering/original/plots/growth-curves/
	python scripts/audano-growth-curves.py -vcf genotype-filtering/augmented/biallelic.filtered.all-samples.vcf.gz -meta igsr_sample_data.tsv -map augmented-graph/map-bub_ids-to-allele_ids.tsv -sample-sheet sample_sheet.tsv -output genotype-filtering/augmented/plots/growth-curves/

boxplots-all:
	python scripts/boxplots-lengths.py -meta igsr_sample_data.tsv -coverages 1k_ont_data_overview-samples_and_sequencing.tsv -names original,original-filtered,augment,augment-filtered -vcfs original-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz,genotype-filtering/original/biallelic.filtered.all-samples.vcf.gz,augmented-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz,genotype-filtering/augmented/biallelic.filtered.all-samples.vcf.gz -output joint-plots/all/boxplot-lengths/
	python scripts/boxplots-counts.py -meta igsr_sample_data.tsv -coverages 1k_ont_data_overview-samples_and_sequencing.tsv -names original,original-filtered,augment,augment-filtered -vcfs original-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz,genotype-filtering/original/biallelic.filtered.all-samples.vcf.gz,augmented-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz,genotype-filtering/augmented/biallelic.filtered.all-samples.vcf.gz -output joint-plots/all/boxplot-counts/

rausch-curves-all:
	python scripts/rausch-curves.py -sample-sheet sample_sheet.tsv -names original,original-filtered,augment,augment-filtered -vcfs original-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz,genotype-filtering/original/biallelic.filtered.all-samples.vcf.gz,augmented-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz,genotype-filtering/augmented/biallelic.filtered.all-samples.vcf.gz -output joint-plots/all/rausch-curves/



# All the scripts for analyzing the phased genotypes

phased-genotypes-analyze: boxplot-phased rausch-curves-phased audano-curves-phased

subset-tag-phased-vcf:
	bcftools view phased-genotypes/phased.biallelic.bcf | awk '($1 ~ /^#/ || $3 != ".")' | bgzip -c > phased-genotypes/phased.biallelic.svs.vcf.gz
	$(CONDA_ACTIVATE) bcftools-latest
	bcftools view phased-genotypes/phased.biallelic.svs.vcf.gz | bcftools +fill-tags | bgzip -c > phased-genotypes/phased.biallelic.svs.tagged.vcf.gz
	$(CONDA_ACTIVATE) base


qc-table-phased:
	python scripts/qc-table-creation-unrelated.py -sample-sheet sample_sheet.tsv -bi-panel augmented-graph/prepared-panel/giggles-ready_biallelic.vcf.gz -multi-panel augmented-graph/prepared-panel/giggles-ready_multiallelic.vcf.gz -bi-callset phased-genotypes/phased.biallelic.svs.tagged.vcf.gz > phased-genotypes/qc.final.tsv

boxplot-phased:
	python scripts/boxplots-lengths.py -meta igsr_sample_data.tsv -coverages 1k_ont_data_overview-samples_and_sequencing.tsv -names phased,augment-filtered-unrel,augment-filtered-all -vcfs phased-genotypes/phased.biallelic.svs.tagged.vcf.gz,genotype-filtering/augmented/biallelic.filtered.unrel-samples.vcf.gz,genotype-filtering/augmented/biallelic.filtered.all-samples.vcf.gz -output joint-plots/phased/boxplot-lengths/
	python scripts/boxplots-counts.py -meta igsr_sample_data.tsv -coverages 1k_ont_data_overview-samples_and_sequencing.tsv -names phased,augment-filtered-unrel,augment-filtered-all -vcfs phased-genotypes/phased.biallelic.svs.tagged.vcf.gz,genotype-filtering/augmented/biallelic.filtered.unrel-samples.vcf.gz,genotype-filtering/augmented/biallelic.filtered.all-samples.vcf.gz -output joint-plots/phased/boxplot-counts/

plot-hwe-af-phased:
	python scripts/hwe-af-simple.py -table phased-genotypes/qc.final.tsv -output phased-genotypes/plots/hwe-af/

rausch-curves-phased:
	python scripts/rausch-curves.py -sample-sheet sample_sheet.tsv -names phased,augment-filtered-unrel -vcfs phased-genotypes/phased.biallelic.svs.tagged.vcf.gz,genotype-filtering/augmented/biallelic.filtered.unrel-samples.vcf.gz -output joint-plots/phased/rausch-curves/

audano-curves-phased:
	python scripts/audano-growth-curves.py -vcf phased-genotypes/phased.biallelic.svs.tagged.vcf.gz -meta igsr_sample_data.tsv -map augmented-graph/map-bub_ids-to-allele_ids.tsv -sample-sheet sample_sheet.tsv -output phased-genotypes/plots/growth-curves/


# script for the mendelian consistency

check-mendelian-consistency:
	bash scripts/mendelian-check.sh
	
plot-mendelian-consistency:
	python scripts/mendelian-consistency-upsetplot.py -children NA19828,HG01258,HG00420,NA19129,NA12877,NA12878 -out augmented-graph/mendelian-consistency/
	python scripts/mendelian-consistency-upsetplot.py -children NA19828,HG01258,HG00420,NA19129,NA12877,NA12878 -out original-graph/mendelian-consistency/
	python scripts/mendelian-consistency-upsetplot.py -children NA19828,HG01258,HG00420,NA19129,NA12877,NA12878 -out genotype-filtering/augmented/mendelian-consistency/
	python scripts/mendelian-consistency-upsetplot.py -children NA19828,HG01258,HG00420,NA19129,NA12877,NA12878 -out genotype-filtering/original/mendelian-consistency/