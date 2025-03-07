.ONESHELL:

# Need to specify bash in order for conda activate to work.
SHELL=/bin/bash
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

# Some metadata files provided along with the repository
UNRELATED_SAMPLES=unrel.samples.tsv
IGSR_METADATA=igsr_sample_data.tsv
SAMPLE_METADATA=sample_sheet.tsv
COVERAGES_OVERVIEW=1k_ont_data_overview-samples_and_sequencing.tsv

### Variables that need to be declared before the makefile is run.

# Main project directory where all the files are present and where to put the joint plots
PROJECT_DIR=results/

# These are the genotypes from the genotyping pipeline
# Biallelic panel made in the rule annotate_panel_giggles_ready of the "hhu-snakemake-pipelines/post-augmentation-genotyping/rules/prepare-vcf.smk" file
# File can be downloaded from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.1/giggles-genotyping/giggles-panel-biallelic.vcf.gz
BI_PANEL=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/giggles-genotyping/augmented-graph/prepared-panel/giggles-ready_biallelic.vcf.gz
# Multiallelic panel made in the rule annotate_panel_giggles_ready of the "hhu-snakemake-pipelines/post-augmentation-genotyping/rules/prepare-vcf.smk" file
# File can be downloaded from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.1/giggles-genotyping/giggles-panel-multiallelic.vcf.gz
MULTI_PANEL=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/giggles-genotyping/augmented-graph/prepared-panel/giggles-ready_multiallelic.vcf.gz

# Biallelic genotypes produced by Giggles on the HPRC_mg in the rule "convert_multisample_vcf_biallelic" of "hhu-snakemake-pipelines/pre-augmentation-giggles-svarp/rules/genotyping.smk" file
# File not part of release. Contact authors for access.
ORIGINAL_GIGGLES_GTS_BA=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/giggles-genotyping/original-graph/multisample-genotypes/multisample-biallelic.vcf.gz
# Multiallelic genotypes produced by Giggles on the HPRC_mg in the rule "merge_vcf_to_multisample" of "hhu-snakemake-pipelines/pre-augmentation-giggles-svarp/rules/genotyping.smk" file
# File not part of release. Contact authors for access.
ORIGINAL_GIGGLES_GTS_MA=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/giggles-genotyping/original-graph/multisample-genotypes/multisample-multiallelic.vcf.gz
# The directory where the above two files are present
ORIGINAL_GIGGLES_GTS_DIR=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/giggles-genotyping/original-graph/multisample-genotypes/

# Biallelic genotypes produced by Giggles on the HPRC_mg_44+966 in the rule "convert_multisample_vcf_biallelic" of "hhu-snakemake-pipelines/post-augmentation-genotyping/rules/genotyping.smk" file
# File can be downloaded from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.1/giggles-genotyping/giggles-genotypes-biallelic.vcf.gz
AUGMENTED_GIGGLES_GTS_BA=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/giggles-genotyping/augmented-graph/multisample-genotypes/multisample-biallelic.vcf.gz
# Multiallelic genotypes produced by Giggles on the HPRC_mg_44+966 in the rule "merge_vcf_to_multisample" of "hhu-snakemake-pipelines/post-augmentation-genotyping/rules/genotyping.smk" file
# File can be downloaded from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.1/giggles-genotyping/giggles-genotypes-multiallelic.vcf.gz
AUGMENTED_GIGGLES_GTS_MA=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/giggles-genotyping/augmented-graph/multisample-genotypes/multisample-multiallelic.vcf.gz
#The directory where the above two files are present
AUGMENTED_GIGGLES_GTS_DIR=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/giggles-genotyping/augmented-graph/multisample-genotypes/

# These are the biallelic genotypes from the post-phasing of the genotypes using ShapeIt5
# These files were produced by Dr Tobias Rausch of EMBL and are not part of the pipeline. Please refer to the Methods section on how to generate them.
# File can be downloaded from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.1/shapeit5-phased-callset/shapeit5-phased-callset_final-vcf.phased.vcf.gz
PHASED_GTS=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/final-vcf.phased.vcf.gz
# The directory of the phased genotypes
PHASED_GTS_DIR=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/


### Derived variables. These are files that will be produced as a result.

ORIGINAL_GIGGLES_GTS_FILTERED=$(ORIGINAL_GIGGLES_GTS_DIR)/biallelic.filtered.all-samples.vcf.gz
AUGMENTED_GIGGLES_GTS_FILTERED=$(AUGMENTED_GIGGLES_GTS_DIR)/biallelic.filtered.all-samples.vcf.gz

ORIGINAL_MAP_ID=$(ORIGINAL_GIGGLES_GTS_DIR)/map-bub_ids-to-allele_ids.tsv
AUGMENTED_MAP_ID=$(AUGMENTED_GIGGLES_GTS_DIR)/map-bub_ids-to-allele_ids.tsv

PHASED_GTS_TAGGED=$(PHASED_GTS_DIR)/phased.gts.tagged.vcf.gz

# QC Tables
PHASED_GTS_QC_TABLE=$(PHASED_GTS_DIR)/phased.qc.tsv
AUGMENTED_GIGGLES_GTS_FILTERED_QC_TABLE=$(AUGMENTED_GIGGLES_GTS_DIR)/filtered.qc.tsv

all: directory-tree map-ids filtering tag-phased-vcf qc-table boxplot rausch-curves hwe-af growth-curves venn-diagram mendelian-consistency get-statistics

# removes all the files and folders made during the execution of the script
remove:
	rm -r $(PROJECT_DIR)/statistics $(PROJECT_DIR)/plots $(PHASED_GTS_DIR)/plots $(AUGMENTED_GIGGLES_GTS_DIR)/plots $(ORIGINAL_GIGGLES_GTS_DIR)/plots $(AUGMENTED_GIGGLES_GTS_DIR)/mendelian-consistency
	rm sample-lists/*.svg
	rm $(ORIGINAL_GIGGLES_GTS_FILTERED)
	rm $(AUGMENTED_GIGGLES_GTS_FILTERED)
	rm $(ORIGINAL_MAP_ID)
	rm $(AUGMENTED_MAP_ID)
	rm $(PHASED_GTS_TAGGED)
	rm $(PHASED_GTS_QC_TABLE)
	rm $(AUGMENTED_GIGGLES_GTS_FILTERED_QC_TABLE)

# creates the inital directory structure where the output files can be stored
directory-tree:
	mkdir -p $(PROJECT_DIR)
	mkdir -p $(PROJECT_DIR)/statistics
	mkdir -p $(PROJECT_DIR)/plots
	mkdir -p $(PHASED_GTS_DIR)/plots
	mkdir -p $(AUGMENTED_GIGGLES_GTS_DIR)/plots
	mkdir -p $(ORIGINAL_GIGGLES_GTS_DIR)/plots
	mkdir -p $(AUGMENTED_GIGGLES_GTS_DIR)/mendelian-consistency

# creating a map between the bubble ids and allele ids
map-ids:
	python scripts/map-bub-ids-to-allele-ids.py -vcf $(ORIGINAL_GIGGLES_GTS_MA) > $(ORIGINAL_MAP_ID)
	python scripts/map-bub-ids-to-allele-ids.py -vcf $(AUGMENTED_GIGGLES_GTS_MA) > $(AUGMENTED_MAP_ID)

# bcftools based filtering
filtering:
	bcftools +fill-tags $(AUGMENTED_GIGGLES_GTS_BA) -- -t all | bcftools filter -i 'INFO/AC >= 1 && INFO/F_MISSING <= 0.05' | bcftools view --min-ac 1 -m2 -M2 | bcftools +setGT - -- -t q -n . -i 'FMT/GQ<10' | bcftools +fill-tags - -- -t all | bcftools filter -O z -o $(AUGMENTED_GIGGLES_GTS_FILTERED) -i 'INFO/AC >= 1 && INFO/F_MISSING <= 0.05' -
	bcftools +fill-tags $(ORIGINAL_GIGGLES_GTS_BA) -- -t all | bcftools filter -i 'INFO/AC >= 1 && INFO/F_MISSING <= 0.05' | bcftools view --min-ac 1 -m2 -M2 | bcftools +setGT - -- -t q -n . -i 'FMT/GQ<10' | bcftools +fill-tags - -- -t all | bcftools filter -O z -o $(ORIGINAL_GIGGLES_GTS_FILTERED) -i 'INFO/AC >= 1 && INFO/F_MISSING <= 0.05' -
	
# tagging the vcf with bcftools
tag-phased-vcf:
	bcftools view $(PHASED_GTS) | awk '($$1 ~ /^#/ || $$3 != ".")' | bcftools +fill-tags | bgzip -c > $(PHASED_GTS_TAGGED)
	
# creating quality check tables out of the vcfs
qc-table:
	python scripts/qc-table-creation-unrelated.py -sample-sheet $(SAMPLE_METADATA) -bi-panel $(BI_PANEL) -multi-panel $(MULTI_PANEL) -bi-callset $(PHASED_GTS_TAGGED) > $(PHASED_GTS_QC_TABLE)
	python scripts/qc-table-creation.py -sample-sheet $(SAMPLE_METADATA) -bi-panel $(BI_PANEL) -multi-panel $(MULTI_PANEL) -bi-callset $(AUGMENTED_GIGGLES_GTS_FILTERED) > $(AUGMENTED_GIGGLES_GTS_FILTERED_QC_TABLE)

# creating boxplot of SV lengths and SV counts
boxplot:
	python scripts/boxplots-lengths.py -meta $(IGSR_METADATA) -coverages $(COVERAGES_OVERVIEW) -names HPRC_mg,HPRC_mg_44+966 -vcfs $(ORIGINAL_GIGGLES_GTS_FILTERED),$(PHASED_GTS_TAGGED) -output $(PROJECT_DIR)/plots
	python scripts/boxplots-counts.py -meta $(IGSR_METADATA) -coverages $(COVERAGES_OVERVIEW) -names HPRC_mg,HPRC_mg_44+966 -vcfs $(ORIGINAL_GIGGLES_GTS_FILTERED),$(PHASED_GTS_TAGGED) -output $(PROJECT_DIR)/plots

# creating log(allele count) vs log(number of variant sites)
rausch-curves:
	python scripts/rausch-curves.py -sample-sheet $(SAMPLE_METADATA) -names HPRC_mg,HPRC_mg_44+966 -vcfs $(ORIGINAL_GIGGLES_GTS_FILTERED),$(PHASED_GTS_TAGGED) -output $(PROJECT_DIR)/plots

# creating HWE (heterozygosity vs allele frequency) and AF (allele frequency in reference panel vs allele frequency in genotyped samples) plots
hwe-af:
	python scripts/hwe-af-simple.py -table $(PHASED_GTS_QC_TABLE) -output $(PHASED_GTS_DIR)/plots/
	python scripts/hwe-af-simple.py -table $(AUGMENTED_GIGGLES_GTS_FILTERED_QC_TABLE) -output $(AUGMENTED_GIGGLES_GTS_DIR)/plots/

# creating SV growth curves with AFR and non-AFR samples kept separate
growth-curves:
	python scripts/audano-growth-curves.py -vcf $(PHASED_GTS_TAGGED) -meta $(IGSR_METADATA) -map $(AUGMENTED_MAP_ID) -sample-sheet $(SAMPLE_METADATA) -output $(PHASED_GTS_DIR)/plots/
	python scripts/audano-growth-curves.py -vcf $(ORIGINAL_GIGGLES_GTS_FILTERED) -meta $(IGSR_METADATA) -map $(ORIGINAL_MAP_ID) -sample-sheet $(SAMPLE_METADATA) -output $(ORIGINAL_GIGGLES_GTS_DIR)/plots/

# creating Venn Diagram of different consortiums to show overlap of data
venn-diagram:
	python scripts/make-sample-venn.py

# calculating number of mendelian inconsistent sites in the six families present in the data
mendelian-consistency:
	bash scripts/mendelian-check.sh $(AUGMENTED_MAP_ID) $(AUGMENTED_GIGGLES_GTS_FILTERED) $(AUGMENTED_GIGGLES_GTS_DIR)/mendelian-consistency
	python scripts/mendelian-consistency-upsetplot.py -children NA19828,HG01258,HG00420,NA19129,NA12877,NA12878 -out $(AUGMENTED_GIGGLES_GTS_DIR)/mendelian-consistency

# some statistics for the paper
get-statistics:
	python scripts/find-map.py -map $(AUGMENTED_MAP_ID) -vcf $(PHASED_GTS) -outdir $(PROJECT_DIR)/statistics/phased-
	python scripts/find-map.py -map $(AUGMENTED_MAP_ID) -vcf $(AUGMENTED_GIGGLES_GTS_FILTERED) -outdir $(PROJECT_DIR)/statistics/unphased-