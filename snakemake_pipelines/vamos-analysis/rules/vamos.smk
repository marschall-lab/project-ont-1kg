# unzip the bed file
rule unzip_vamos_sites_list:
    input:
        'resources/vamos-sites-list.bed.gz'
    output:
        temp('resources/vamos-sites-list.bed')
    conda:
        '../envs/vamos.yml'
    shell:
        'gzip -dk {input}'

# run vamos on the vienna data
rule vamos_vienna:
    input:
        alignment=config['path_to_cram']+'{sample}.hg38.cram',
        sites='resources/vamos-sites-list.bed'
    output:
        vcf='results/vamos-results-vienna/{sample}.vcf'
    log:
        'results/vamos-results-vienna/{sample}.log'
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=5000
    shell:
        'vamos --read -b {input.alignment} -r {input.sites} -s {wildcards.sample} -o {output} > {log}'