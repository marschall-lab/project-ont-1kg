# 1000 Genomes ONT Sequencing

## SV Annotation by Graph Augmentation (SAGA) framework
![overview_figure](https://github.com/marschall-lab/project-ont-1kg/blob/main/figures/SAGA-framework.png)

The repository hosts the snakemake pipelines and analysis scripts used for analysing the 1,019 samples sequenced with ONT long-reads.

The following snakemake pipelines are hosted in this repository:

- Haplotagging of Aligned Reads under `haplotagging`: The pipeline tags the aligned reads of this study as originating from haplotype 1 or 2 using `whatshap`.
- Phasing Experiments under `phasing`: The pipeline phases the NYGC raw genotypes using the aligned reads with `whatshap` and we perform QC using the NYGC statistical phased VCF.
- Running Giggles and SVarp on the HPRC_mg graph under `pre-augmentation-giggles-svarp`: The pipeline pre-processes and does SV discovery with `SVarp` and SV genotyping with `giggles` using HPRC_mg.
- Running Giggles on the HPRC_mg_44+966 graph under `post-augmentation-genotyping`: The pipeline pre-processes and does SV genotyping with `giggles` using HPRC_mg_44+966.
- Running QCs based on coverage and read N50 stratification under `genotype-sample-subset-analysis`: The pipeline bins the samples based on coverage and read N50 and runs QC to check the effects of these two variables on genotyping and SV discovery.
- Investigating recovery of added SVs during genotyping under `sv-recovery`: The pipeline investigates the SVs that were introduced into the graph to create HPRC_mg_44+966 and whether they were retrieved during genotyping.
- Annotating the alleles in HPRC_mg_44+966 based on the ancestral allele found in Chimpanzee under `chimpanzee-analysis`: The pipeline identifies ancestral allele in the chimpanzee using alignment of the chimpanzee assembly to HPRC_mg_44+966 and tags the VCFs produced from `giggles`.
- VNTR genotyping of the cohort using vamos under `vamos-analysis`: The pipeline runs the VNTR calling using `vamos` and subsequent QC with the HGSVC3 assemblies.

For information about the other scripts, please refer to the "Code Availability" section of the preprint.

Link to preprint: https://www.biorxiv.org/content/10.1101/2024.04.18.590093v1
