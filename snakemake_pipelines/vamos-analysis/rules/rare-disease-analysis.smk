
# plot histogram
rule disease_vntr_histogram:
    input:
        ont=expand('results/vamos-t2t/{sample}.stats', sample=cram_sample_list),
        hgsvc=expand('results/hgsvc3-comparison/vntr-calls/{sample}.stats', sample=hgsvc_sample_list_all)
    output:
        'results/rare-disease-vntr-analysis/abca7-histogram.svg',
        'results/rare-disease-vntr-analysis/plin4-histogram.svg'
    params:
        outdir='results/rare-disease-vntr-analysis/'
        hgsvc=','.join(list(expand('results/hgsvc3-comparison/vntr-calls/{sample}.stats', sample=hgsvc_sample_list_all))),
        ont=','.join(list(expand('results/vamos-t2t/{sample}.vcf', sample=cram_sample_list)))
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=5,
        runtime_min=0,
        mem_total_mb=10*1024
    shell:
        'python scripts/plot-vntr-disease-histogram.py -hgsvc {params.hgsvc} -ont {params.ont} -output {params.outdir}'