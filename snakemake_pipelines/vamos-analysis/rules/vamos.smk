## Running VAMOS on hg38 for comparison against Miller vcf

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

rule sort_bed_file:
    input:
        'resources/vamos-sites-list.bed'
    output:
        temp('resources/vamos-sites-list.sorted.bed')
    shell:
        'sort -k 1,1 -k2,2n {input} > {output}'

# convert crams to bams
rule cram_to_bam_hg38:
    input:
        reads=config['path_to_hg38_cram']+'{sample}.hg38.cram',
        ref=config['reference_directory']+'1KG_ONT_VIENNA_hg38.fa'
    output:
        bam=temp('results/temp/{sample}.hg38.bam'),
        ind=temp('results/temp/{sample}.hg38.bam.bai')
    params:
        ref_dir=config['reference_directory']
    conda:
        '../envs/vamos.yml'
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
        alignment='results/temp/{sample}.hg38.bam',
        alignment_index='results/temp/{sample}.hg38.bam.bai',
        sites='resources/vamos-sites-list.sorted.bed'
    output:
        vcf='results/miller-comparison/vamos-results-vienna/{sample}.vcf'
    log:
        'results/miller-comparison/vamos-results-vienna/{sample}.log'
    threads: 8
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=24,
        runtime_min=20,
        mem_total_mb=500*1024
    shell:
        'vamos --read -b {input.alignment} -r {input.sites} -s {wildcards.sample} -o {output} -t {threads} > {log}'

# process sites list to make it compatible with bedtools
rule process_sites_list:
    input:
        'resources/vamos-sites-list.sorted.bed'
    output:
        temp('resources/vamos-sites-list.sorted.first3columns.bed')
    shell:
        'cat {input} | cut -f 1-3 > {output}'

# extract VNTR sequence from the GRCh38 reference
rule extract_reference_sequence:
    input:
        ref=config['reference_directory']+'1KG_ONT_VIENNA_hg38.fa',
        sites='resources/vamos-sites-list.sorted.first3columns.bed'
    output:
        'results/miller-comparison/reference-vntrs.bed'
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=5*1024
    shell:
        'bedtools getfasta -bedOut -fi {input.ref} -bed {input.sites} > {output}'
        

## Running Vamos cohort-wide on T2T

# processing T2T sites list
# the VAMOS list was obtained from https://zenodo.org/records/13263615
rule process_t2t_sites_list:
    input:
        'resources/vamos.effMotifs-0.1.T2T-CHM13.tsv.gz'
    output:
        'results/temp/vamos.T2T.processed.tsv'
    shell:
        'zcat {input} | grep VNTR | sort -k 1,1 -k2,2n > {output}'

# convert crams to bams
rule cram_to_bam_t2t:
    input:
        reads=config['path_to_t2t_cram']+'{sample}.t2t.cram',
        ref=config['reference_directory']+'1KG_ONT_VIENNA_t2t.fa'
    output:
        bam=temp('results/temp/{sample}.t2t.bam'),
        ind=temp('results/temp/{sample}.t2t.bam.bai')
    params:
        ref_dir=config['reference_directory']
    conda:
        '../envs/vamos.yml'
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

# run vamos on the vienna data cohort-wide
rule vamos_t2t:
    input:
        alignment='results/temp/{sample}.t2t.bam',
        alignment_index='results/temp/{sample}.t2t.bam.bai',
        sites='results/temp/vamos.T2T.processed.tsv'
    output:
        vcf='results/vamos-t2t/{sample}.vcf'
    log:
        'results/vamos-t2t/{sample}.log'
    threads: 8
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=12,
        runtime_min=20,
        mem_total_mb=50*1024
    shell:
        'vamos --read -b {input.alignment} -r {input.sites} -s {wildcards.sample} -o {output} -t {threads} > {log}'

# getting stats from the vntr vcf
rule vamos_t2t_stats:
    input:
        'results/vamos-t2t/{sample}.vcf'
    output:
        'results/vamos-t2t/{sample}.stats'
    shell:
        'python scripts/get-vntr-stats.py -vcf {input} > {output}'

# create vntr summary table
rule vamos_t2t_summary:
    input:
        vcf=expand('results/vamos-t2t/{sample}.stats', sample=t2t_samples),
        sites='results/temp/vamos.T2T.processed.tsv'
    output:
        'results/vamos-t2t-summary.bed'
    params:
        ','.join(list(expand('results/vamos-t2t/{sample}.stats', sample=t2t_samples)))
    conda:
        '../envs/comparison.yml'
    resources:
        runtime_hrs=5,
        runtime_min=0,
        mem_total_mb=96000
    shell:
        'python scripts/create-vntr-table.py -stats {params} -sites {input.sites} > {output}'

# combine sample vcfs into one multisample vcf
rule vamos_t2t_vcf_combine:
    input:
        vcf=expand('results/vamos-t2t/{sample}.vcf', sample=t2t_samples)
    output:
        'results/vamos-multisample.vcf'
    params:
        ','.join(list(expand('results/vamos-t2t/{sample}.vcf', sample=t2t_samples)))
    shell:
        'python scripts/combine-vamos-vcf.py -vcfs {params} > {output}'