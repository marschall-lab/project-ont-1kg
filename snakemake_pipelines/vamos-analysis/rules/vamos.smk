# unzip the bed file
rule unzip_vamos_sites_list:
    input:
        'resources/vamos-sites-list.bed.gz'
    output:
        temp('resources/vamos-sites-list.bed')
    shell:
        'gzip -dk {input}'

# run vamos on the vienna data
rule vamos_vienna:
    input:
        alignments=config['path_to_cram']+'/{sample}.hg38.cram',
        sites='resources/vamos-sites-list.bed'
    output:
        vcf='results/vamos-results-vienna/{sample}.vcf'
    log:
        'results/vamos-results-vienna/{sample}.log'
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=0,
        runtime_min=20,
        mem_total_mb=256
    shell:
        'vamos --read -b {input.alignment} -r {input.sites} -s {wildcard.sample} -o {output} > {log}'