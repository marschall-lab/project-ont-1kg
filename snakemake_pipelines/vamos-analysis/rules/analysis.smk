## Adding rules for creating BED files for each sample and show which VNTR is non-reference.

# align chm13 assembly to itself
rule align_chm13_to_itself:
    input:
        assembly=config['reference_directory']+'1KG_ONT_VIENNA_t2t.fa',
    output:
        tmp_sam=temp('results/per-sample-bed/chm13-prep/tmp.sam'),
        alignment='results/per-sample-bed/chm13-prep/chm13_mapped-to-chm13.bam'
    conda:
        '../envs/prepare-bed.yml'
    resources:
        runtime_hrs=12,
        runtime_min=20,
        mem_total_mb=50*1024
    shell:
        '''
        minimap2 -ax asm5 {input.assembly} {input.assembly} > {output.tmp_sam}
        samtools view -bS {output.tmp_sam} > {output.alignment}
        '''

# running vamos in contig mode
rule call_chm13_vntr:
    input:
        alignment='results/per-sample-bed/chm13-prep/chm13_mapped-to-chm13.bam',
        vntr_sites='results/temp/vamos.T2T.processed.tsv'
    output:
        'results/per-sample-bed/chm13-prep/chm13-vntr.vcf'
    log:
        'results/per-sample-bed/chm13-prep/chm13-vntr.log'
    threads: 8
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=12,
        runtime_min=20,
        mem_total_mb=50*1024
    shell:
        'vamos --contig -b {input.alignment} -r {input.vntr_sites} -s chm13 -o {output} -t {threads} > {log}'

# get stats from the ref vntrs from vamos
rule chm13_vntr_stats:
    input:
        vcf='results/per-sample-bed/chm13-prep/chm13-vntr.vcf'
    output:
        stats='results/per-sample-bed/chm13-prep/chm13-vntr.stats',
        plot='results/per-sample-bed/chm13-prep/chm13-vntr-length-diff.png'
    conda:
        '../envs/vamos.yml'
    shell:
        'python scripts/get-ref-vntr-stats.py -vcf {input.vcf} -plot {output.plot} > {output.stats}'    


# creating BED file for sample with its non-reference records.
rule create_samplewise_vntr_bed:
    input:
        sample_vntr='results/vamos-t2t/{sample}.vcf',
        ref_vntr='results/per-sample-bed/chm13-prep/chm13-vntr.vcf'
    output:
        bed='results/per-sample-bed/bed/{sample}.bed',
        stats='results/per-sample-bed/bed/{sample}.stats'
    resources:
        runtime_hrs=1,
        runtime_min=0,
        mem_total_mb=1024
    shell:
        'python scripts/prepare-bed.py -sample {input.sample_vntr} -ref {input.ref_vntr} 1> {output.bed} 2> {output.stats}'


# extract VNTR sequence from the CHM13-T2T reference
rule extract_reference_sequence:
    input:
        ref=config['reference_directory']+'1KG_ONT_VIENNA_t2t.fa',
        sites='results/temp/vamos.T2T.processed.tsv'
    output:
        sites=temp('results/temp/vamos-sites-list.sorted.first3columns.bed'),
        bed='results/per-sample-bed/reference-vntrs-t2t-seq.bed'
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=5*1024
    shell:
        '''
        cat {input.sites} | cut -f 1-3 > {output.sites}
        bedtools getfasta -bedOut -fi {input.ref} -bed {input.sites} > {output.bed}
        '''

# creating VCF file for sample.
rule create_samplewise_vntr_vcf:
    input:
        sample_vntr='results/vamos-t2t/{sample}.vcf',
        ref_vntr='results/per-sample-bed/chm13-prep/chm13-vntr.vcf',
        ref_seq='results/per-sample-bed/reference-vntrs-t2t-seq.bed'
    output:
        'results/per-sample-bed/vcf/{sample}.vcf.gz'
    resources:
        runtime_hrs=1,
        runtime_min=0,
        mem_total_mb=1024
    shell:
        'python scripts/prepare-vntr-vcf.py -sample {input.sample_vntr} -ref {input.ref_vntr} -refseq {input.ref_seq} | bgzip -c > {output}'

# curating the stats
rule curate_samplewise_vntr_stats:
    input:
        stats=expand('results/per-sample-bed/bed/{sample}.stats', sample=cram_sample_list)
    output:
        'results/per-sample-bed/bed-stats.tsv'
    params:
        ','.join(list(expand('results/per-sample-bed/bed/{sample}.stats', sample=cram_sample_list)))
    shell:
        'python scripts/curate-bed-stats.py -stats {params} > {output}'