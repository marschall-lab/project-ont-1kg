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

# convert crams to bams
rule cram_to_bam:
    input:
        reads=config['path_to_cram']+'/{sample}.hg38.cram',
        ref=config['reference_directory']+'/1KG_ONT_VIENNA_hg38.fa'
    output:
        bam=temp('results/temp/{sample}.bam'),
        ind=temp('results/temp/{sample}.bam.bai')
    params:
        ref_dir=config['reference_directory']
    conda:
        'envs/vamos.yml'
    resources:
        runtime_hrs=12,
        runtime_min=0,
        mem_total_mb=5000
    shell:
        '''
        seq_cache_populate.pl -root {params.ref_dir} {input.ref}
        export REF_PATH={params.ref_dir}/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
        export REF_CACHE={params.ref_dir}/%2s/%2s/%s
        samtools view -b -T {input.ref} -o {output.bam} {input.reads}
        samtools index -b {output.bam}
        '''

# run vamos on the vienna data
rule vamos_vienna:
    input:
        alignment='results/temp/{sample}.bam',
        sites='resources/vamos-sites-list.bed'
    output:
        vcf='results/vamos-results-vienna/{sample}.vcf'
    log:
        'results/vamos-results-vienna/{sample}.log'
    threads: 24
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=24,
        runtime_min=20,
        mem_total_mb=5000
    shell:
        'vamos --read -b {input.alignment} -r {input.sites} -s {wildcards.sample} -o {output} -t {threads} > {log}'