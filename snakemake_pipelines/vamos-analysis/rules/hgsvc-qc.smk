## Running QC using the HGSVC assemblies

# running VNTR calling on assemblies.
rule hgsvc_vntr_call:
    input:
        bam=config['path_to_hgsvc_alignments']+'{sample}.vrk-ps-sseq.asm-{haplotype}.t2tv2.sort.bam',
        sites='results/temp/vamos.T2T.processed.tsv'
    output:
        'results/hgsvc3-comparison/vntr-calls/{sample}-{haplotype}.vcf'
    log:
        'results/hgsvc3-comparison/vntr-calls/{sample}-{haplotype}.log'
    threads: 8
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=12,
        runtime_min=20,
        mem_total_mb=50*1024
    shell:
        'vamos --contig -b {input.bam} -r {input.tsv} -s {wildcards.sample} -o {output} -t {threads}'