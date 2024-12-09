.ONESHELL:

# Need to specify bash in order for conda activate to work.
SHELL=/bin/bash
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

# Main project directory where all the files are present and where to put the results plots
PROJECT_DIR=results

# These are the genotypes from the genotyping pipeline
# Biallelic panel made in the rule annotate_panel_giggles_ready of the "hhu-snakemake-pipelines/post-augmentation-genotyping/rules/prepare-vcf.smk" file
BI_PANEL=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/giggles-genotyping/augmented-graph/prepared-panel/giggles-ready_biallelic.vcf.gz
# Multiallelic panel made in the rule annotate_panel_giggles_ready of the "hhu-snakemake-pipelines/post-augmentation-genotyping/rules/prepare-vcf.smk" file
MULTI_PANEL=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/giggles-genotyping/augmented-graph/prepared-panel/giggles-ready_multiallelic.vcf.gz

# Multiallelic genotypes produced by Giggles on the HPRC_mg_44+966 in the rule "merge_vcf_to_multisample" of "hhu-snakemake-pipelines/post-augmentation-genotyping/rules/genotyping.smk" file
# The above pipeline produced file has been processed using the analysis_scripts/genotyping Makefile
AUGMENTED_GIGGLES_GTS_FILTERED=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/giggles-genotyping/augmented-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz

# path to the SVAN insertion annotations provided by Berni
SVAN_INS=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/SVAN/29_11_2024_resubmission/final-vcf.unphased.INS.SVAN_1.3.fixed.vcf

# path to the SVAN deletion annotations provided by Berni
SVAN_DEL=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/SVAN/29_11_2024_resubmission/final-vcf.unphased.DEL.SVAN_1.3.fixed.vcf

# path to the SVAN complex insertion annotations provided by Berni
SVAN_COMPLEX_INS=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/SVAN/29_11_2024_resubmission/final-vcf.unphased.COMPLEX-INS.SVAN_1.3.fixed.vcf

# path to the SVAN complex deletion annotations provided by Berni
SVAN_COMPLEX_DEL=/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/SVAN/29_11_2024_resubmission/final-vcf.unphased.COMPLEX-DEL.SVAN_1.3.fixed.vcf

# path to the list of VNTR sites given to vamos for the pipeline in snakemake_pipelines/vamos-analysis.
# the output of the rule process_t2t_sites_list in rules/vamos.smk
VAMOS_SITES=/gpfs/project/projects/medbioinf/users/spani/scripts/1000g-ont/ont-1kg/snakemake_pipelines/vamos-analysis/results/temp/vamos.T2T.processed.tsv

all: directory-tree gts extract-categories rausch-one qc vntr-stratification intersect-with-vamos

# removes all the files and folders made during the execution of the script
remove:
	rm -r $(PROJECT_DIR)
	
# creates the inital directory structure where the output files can be stored
directory-tree:
	mkdir -p $(PROJECT_DIR)
	mkdir -p $(PROJECT_DIR)/vcfs
	mkdir -p $(PROJECT_DIR)/tables
	mkdir -p $(PROJECT_DIR)/statistics

gts:
	bcftools sort $(SVAN_INS) | bgzip -c > $(PROJECT_DIR)/vcfs/svan-ins-sorted.vcf.gz
	bcftools sort $(SVAN_DEL) | bgzip -c > $(PROJECT_DIR)/vcfs/svan-del-sorted.vcf.gz
	bcftools sort $(SVAN_COMPLEX_INS) | bgzip -c > $(PROJECT_DIR)/vcfs/svan-complex-ins-sorted.vcf.gz
	bcftools sort $(SVAN_COMPLEX_DEL) | bgzip -c > $(PROJECT_DIR)/vcfs/svan-complex-del-sorted.vcf.gz

	tabix -p vcf $(PROJECT_DIR)/vcfs/svan-ins-sorted.vcf.gz
	tabix -p vcf $(PROJECT_DIR)/vcfs/svan-del-sorted.vcf.gz
	tabix -p vcf $(PROJECT_DIR)/vcfs/svan-complex-ins-sorted.vcf.gz
	tabix -p vcf $(PROJECT_DIR)/vcfs/svan-complex-del-sorted.vcf.gz

	bcftools annotate -c INFO -a $(PROJECT_DIR)/vcfs/svan-ins-sorted.vcf.gz $(AUGMENTED_GIGGLES_GTS_FILTERED) | bcftools sort | bgzip -c > $(PROJECT_DIR)/vcfs/INS.tmp.vcf.gz
	bcftools annotate -c INFO -a $(PROJECT_DIR)/vcfs/svan-del-sorted.vcf.gz $(AUGMENTED_GIGGLES_GTS_FILTERED) | bcftools sort | bgzip -c > $(PROJECT_DIR)/vcfs/DEL.tmp.vcf.gz
	bcftools annotate -c INFO -a $(PROJECT_DIR)/vcfs/svan-complex-ins-sorted.vcf.gz $(AUGMENTED_GIGGLES_GTS_FILTERED) | bcftools sort | bgzip -c > $(PROJECT_DIR)/vcfs/COMPLEX-INS.tmp.vcf.gz
	bcftools annotate -c INFO -a $(PROJECT_DIR)/vcfs/svan-complex-del-sorted.vcf.gz $(AUGMENTED_GIGGLES_GTS_FILTERED) | bcftools sort | bgzip -c > $(PROJECT_DIR)/vcfs/COMPLEX-DEL.tmp.vcf.gz
	
	tabix -p vcf $(PROJECT_DIR)/vcfs/INS.tmp.vcf.gz
	tabix -p vcf $(PROJECT_DIR)/vcfs/DEL.tmp.vcf.gz
	tabix -p vcf $(PROJECT_DIR)/vcfs/COMPLEX-INS.tmp.vcf.gz
	tabix -p vcf $(PROJECT_DIR)/vcfs/COMPLEX-DEL.tmp.vcf.gz

	bcftools concat -a $(PROJECT_DIR)/vcfs/INS.tmp.vcf.gz $(PROJECT_DIR)/vcfs/COMPLEX-INS.tmp.vcf.gz > $(PROJECT_DIR)/vcfs/INS.vcf
	bcftools concat -a $(PROJECT_DIR)/vcfs/DEL.tmp.vcf.gz $(PROJECT_DIR)/vcfs/COMPLEX-DEL.tmp.vcf.gz > $(PROJECT_DIR)/vcfs/DEL.vcf

	rm $(PROJECT_DIR)/vcfs/INS.tmp.vcf.gz $(PROJECT_DIR)/vcfs/DEL.tmp.vcf.gz $(PROJECT_DIR)/vcfs/COMPLEX-INS.tmp.vcf.gz $(PROJECT_DIR)/vcfs/COMPLEX-DEL.tmp.vcf.gz

extract-categories:
	bcftools view -i "ITYPE_N='VNTR'" $(PROJECT_DIR)/vcfs/INS.vcf > $(PROJECT_DIR)/vcfs/VNTR.expansion.gts.vcf
	bcftools view -i "DTYPE_N='VNTR'" $(PROJECT_DIR)/vcfs/DEL.vcf > $(PROJECT_DIR)/vcfs/VNTR.contraction.gts.vcf
	bcftools view -i "ITYPE_N='COMPLEX_DUP' | ITYPE_N='DUP' | ITYPE_N='INV_DUP' | ITYPE_N='DUP_INTERSPERSED'" $(PROJECT_DIR)/vcfs/INS.vcf > $(PROJECT_DIR)/vcfs/DUP.gts.vcf
	bcftools view -i "ITYPE_N='PSD'" $(PROJECT_DIR)/vcfs/INS.vcf > $(PROJECT_DIR)/vcfs/PSD.gts.vcf
	bcftools view -i "ITYPE_N='NUMT'" $(PROJECT_DIR)/vcfs/INS.vcf > $(PROJECT_DIR)/vcfs/NUMT.gts.vcf
	bcftools view -i "ITYPE_N='solo' | ITYPE_N='partnered' | ITYPE_N='orphan'" $(PROJECT_DIR)/vcfs/INS.vcf > $(PROJECT_DIR)/vcfs/MEI.non-ref.gts.vcf
	bcftools view -i "DTYPE_N='solo' | DTYPE_N='partnered' | DTYPE_N='orphan'" $(PROJECT_DIR)/vcfs/DEL.vcf > $(PROJECT_DIR)/vcfs/MEI.ref.gts.vcf

rausch-separate:
	for path in $(PROJECT_DIR)/vcfs/*.gts.vcf; do\
		name=$$(echo $$(path) | cut -d'/' -f 2);\
		echo $${name};\
		python scripts/rausch-curves.py -vcfs $$(path) -names $$(name) -sample-sheet sample_sheet.tsv -output $(PROJECT_DIR)/$$(name);\
	done
	
rausch-one:
	delim="";\
	paths="";\
	names="";\
	for path in $(PROJECT_DIR)/vcfs/*.gts.vcf; do\
		name=$(echo $$path | cut -d'/' -f 2);\
		paths="$(paths)$(delim)$(path)";\
		names="$(names)$(delim)$(name)";\
		delim=",";\
	done;\
	python scripts/rausch-curves.py -vcfs $(paths) -names $(names) -sample-sheet sample_sheet.tsv -output $(PROJECT_DIR)/all;\
	
qc:
	for path in $(PROJECT_DIR)/vcfs/*.gts.vcf; do\
		name=$(echo $$path | cut -d'/' -f 2);\
		echo ${name};\
		python scripts/qc-table-creation.py -sample-sheet sample_sheet.tsv -bi-panel ${BI_PANEL} -multi-panel ${MULTI_PANEL} -bi-callset $(path) > $(PROJECT_DIR)/tables/$$name.tsv;\
	done

# breaking down VNTRs from SVAN as multiallelic or biallelic
vntr-stratification:
	python scripts/map-bub-ids-to-allele-ids.py -vcf $(MULTI_PANEL) > $(PROJECT_DIR)/map-bub_ids-to-allele_ids.tsv
	python scripts/find-vntr-map.py -con $(PROJECT_DIR)/vcfs/VNTR.contraction.gts.vcf -exp $(PROJECT_DIR)/vcfs/VNTR.expansion.gts.vcf -map $(PROJECT_DIR)/map-bub_ids-to-allele_ids.tsv > $(PROJECT_DIR)/statistics/vntr-map.stats

intersect-with-vamos:
	bcftools view -i "ITYPE_N='DUP'" $(PROJECT_DIR)/vcfs/INS.vcf > $(PROJECT_DIR)/statistics/DUP.gts.vcf
	bedtools intersect -a $(VAMOS_SITES) -b $(PROJECT_DIR)/statistics/DUP.gts.vcf > $(PROJECT_DIR)/statistics/dup-instersection.bed
	bedtools intersect -a $(VAMOS_SITES) -b $(PROJECT_DIR)/vcfs/VNTR.expansion.gts.vcf > $(PROJECT_DIR)/statistics/vntr-exp-instersection.bed
	bedtools intersect -a $(VAMOS_SITES) -b $(PROJECT_DIR)/vcfs/VNTR.contraction.gts.vcf > $(PROJECT_DIR)/statistics/vntr-con-instersection.bed