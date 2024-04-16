# Details of the Pipeline

- This snakemake pipeline is to genotype and call variants using the longread dataset and the HPRC_mg graph.
- Genotypes are inferred with Giggles and SVs are called with SVarp.
- The SV calls of SVarp from this pipeline are then input to the graph augmentation pipeline to create the HPRC_mg_44+966 graph.
- The Giggles genotypes here act as a benchmark for comparison to the genotypes on HPRC_mg_44+966.