.ONESHELL:

# Need to specify bash in order for conda activate to work.
SHELL=/bin/bash
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

# Main project directory where all the files are present and where to put the results plots
PROJECT_DIR="results/"

# These are the genotypes from the genotyping pipeline
# Biallelic panel made in the rule annotate_panel_giggles_ready of the "hhu-snakemake-pipelines/post-augmentation-genotyping/rules/prepare-vcf.smk" file
BI_PANEL="/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/giggles-genotyping/augmented-graph/prepared-panel/giggles-ready_biallelic.vcf.gz"
# Multiallelic panel made in the rule annotate_panel_giggles_ready of the "hhu-snakemake-pipelines/post-augmentation-genotyping/rules/prepare-vcf.smk" file
MULTI_PANEL="/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/giggles-genotyping/augmented-graph/prepared-panel/giggles-ready_multiallelic.vcf.gz"

# Multiallelic genotypes produced by Giggles on the HPRC_mg_44+966 in the rule "merge_vcf_to_multisample" of "hhu-snakemake-pipelines/post-augmentation-genotyping/rules/genotyping.smk" file
# The above pipeline produced file has been processed using the analysis_scripts/genotyping Makefile
AUGMENTED_GIGGLES_GTS_FILTERED="/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/giggles-genotyping/augmented-graph/multisample-genotypes/biallelic.filtered.all-samples.vcf.gz"

# path to the SVAN insertion annotations provided by Berni
SVAN_INS="/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/SVAN/29_11_2024_resubmission/final-vcf.unphased.INS.SVAN_1.3.vcf"

# path to the SVAN deletion annotations provided by Berni
SVAN_DEL="/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/SVAN/29_11_2024_resubmission/final-vcf.unphased.DEL.SVAN_1.3.vcf"

# path to the SVAN complex insertion annotations provided by Berni
SVAN_COMPLEX_INS="/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/SVAN/29_11_2024_resubmission/final-vcf.unphased.COMPLEX-INS.SVAN_1.3.vcf"

# path to the SVAN complex deletion annotations provided by Berni
SVAN_COMPLEX_DEL="/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/SVAN/29_11_2024_resubmission/final-vcf.unphased.COMPLEX-DEL.SVAN_1.3.vcf"

all: directory-tree gts extract-categories rausch-one qc

# removes all the files and folders made during the execution of the script
remove:
    rm -r $(PROJECT_DIR)
    
# creates the inital directory structure where the output files can be stored
directory-tree:
    mkdir $(PROJECT_DIR)
    mkdir $(PROJECT_DIR)/vcfs
    mkdir $(PROJECT_DIR)/tables
    
gts:
    bcftools annotate -c INFO -a $(SVAN_INS) $(AUGMENTED_GIGGLES_GTS_FILTERED) | bcftools sort > $(PROJECT_DIR)/vcfs/INS.tmp.vcf
    bcftools annotate -c INFO -a $(SVAN_DEL) $(AUGMENTED_GIGGLES_GTS_FILTERED) | bcftools sort > $(PROJECT_DIR)/vcfs/DEL.tmp.vcf
    bcftools annotate -c INFO -a $(SVAN_COMPLEX_INS) $(AUGMENTED_GIGGLES_GTS_FILTERED) | bcftools sort > $(PROJECT_DIR)/vcfs/COMPLEX-INS.tmp.vcf
    bcftools annotate -c INFO -a $(SVAN_COMPLEX_DEL) $(AUGMENTED_GIGGLES_GTS_FILTERED) | bcftools sort > $(PROJECT_DIR)/vcfs/COMPLEX-DEL.tmp.vcf
    bcftools concat $(PROJECT_DIR)/vcfs/INS.tmp.vcf $(PROJECT_DIR)/vcfs/COMPLEX-INS.tmp.vcf > $(PROJECT_DIR)/vcfs/INS.vcf
    bcftools concat $(PROJECT_DIR)/vcfs/DEL.tmp.vcf $(PROJECT_DIR)/vcfs/COMPLEX-DEL.tmp.vcf > $(PROJECT_DIR)/vcfs/DEL.vcf

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
        name=$$(echo $$path | cut -d'/' -f 2);\
        echo $${name};\
        python scripts/rausch-curves.py -vcfs $$path -names $$name -sample-sheet sample_sheet.tsv -output $(PROJECT_DIR)/$$name;\
    done
    
rausch-one:
    delim="";\
    paths="";\
    names="";\
    for path in vcfs/*.gts.vcf; do\
        name=$$(echo $$path | cut -d'/' -f 2);\
        paths="$$paths$$delim$$path";\
        names="$$names$$delim$$name";\
        delim=",";\
    done;\
    python scripts/rausch-curves.py -vcfs $$paths -names $$names -sample-sheet sample_sheet.tsv -output $(PROJECT_DIR)/all;\
    
qc:
    for path in vcfs/*.gts.vcf; do\
        name=$$(echo $$path | cut -d'/' -f 2);\
        echo $${name};\
        python scripts/qc-table-creation.py -sample-sheet sample_sheet.tsv -bi-panel ${BI_PANEL} -multi-panel ${MULTI_PANEL} -bi-callset $$path > $(PROJECT_DIR)/tables/$$name.tsv;\
    done
