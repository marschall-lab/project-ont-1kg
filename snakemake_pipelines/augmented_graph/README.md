# Description of the Pipeline

- Uses the pseudohaplotypes from the SAGA pipeline, HPRC assemblies, and the rGFA file (which has the pseudohaplotypes included) to create genotype calls from Giggles.
- Creates the panel VCF by aligning assemblies (and pseudohaplotypes) back to the rGFA and finding the alleles from each bubble.
- Aligns the long read DNA sequences to the rGFA using minigraph and then sorts it with gaftools.
- Genotypes the samples from the callset using the panel VCF as the set of genotypable variant positions using Giggles.
- Post-processes some basic graphs and statistics from the genotypes (most of the post-genotyping analysis is shown in the folder `analysis_scripts`)