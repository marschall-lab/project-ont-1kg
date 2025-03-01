# sort miller vcf
rule sort_miller_vcf:
    input:
        config['path_to_miller_vcf']
    output:
        temp('results/miller-comparison/analysis/vamos-miller.sorted.vcf')
    shell:
        "zcat {input} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}"

# create TSV of miller vcf
rule parse_miller_vcf:
    input:
        'results/miller-comparison/analysis/vamos-miller.sorted.vcf'
    output:
        'results/miller-comparison/analysis/miller-tsv/{sample}-count.tsv',
        'results/miller-comparison/analysis/miller-tsv/{sample}-sequence.tsv'
    params:
        out_prefix='results/miller-comparison/analysis/miller-tsv/{sample}'
    shell:
        'python scripts/vamos-parse-miller.py -vcf {input} -sample {wildcards.sample} -output {params.out_prefix}'

# create TSV of vienna vcf
rule parse_vienna_vcf:
    input:
        'results/miller-comparison/vamos-results-vienna/{sample}.vcf'
    output:
        'results/miller-comparison/analysis/vienna-tsv/{sample}-count.tsv',
        'results/miller-comparison/analysis/vienna-tsv/{sample}-sequence.tsv'
    params:
        out_prefix='results/miller-comparison/analysis/vienna-tsv/{sample}'
    shell:
        'python scripts/vamos-parse-vienna.py -vcf {input} -sample {wildcards.sample} -output {params.out_prefix}'

# density plot to compare the count of repeating units in both vcfs
rule ru_densityplot:
    input:
        miller='results/miller-comparison/analysis/miller-tsv/{sample}-count.tsv',
        vienna='results/miller-comparison/analysis/vienna-tsv/{sample}-count.tsv'
    output:
        'results/miller-comparison/plots/ru-density/{sample}.png'
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=1,
        runtime_min=0,
        mem_total_mb=5*1024
    shell:
        'python scripts/plot-ru-count-density.py -miller {input.miller} -vienna {input.vienna} -output {output}'

# density plot to compare the count of repeating units in wrong sample combinations
rule ru_densityplot_cross:
    input:
        miller='results/miller-comparison/analysis/miller-tsv/{sample1}-count.tsv',
        vienna='results/miller-comparison/analysis/vienna-tsv/{sample2}-count.tsv'
    output:
        'results/miller-comparison/plots/ru-density/{sample1}_{sample2}.png'
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=1,
        runtime_min=0,
        mem_total_mb=5*1024
    shell:
        'python scripts/plot-ru-count-density.py -miller {input.miller} -vienna {input.vienna} -output {output}'

# histogram plot to compare the allele distance between VNTRs
rule distance_histogram:
    input:
        miller='results/miller-comparison/analysis/miller-tsv/{sample}-sequence.tsv',
        vienna='results/miller-comparison/analysis/vienna-tsv/{sample}-sequence.tsv',
        reference='results/miller-comparison/reference-vntrs.bed'
    output:
        'results/miller-comparison/plots/distance-histogram/{sample}-{af}-subset100.png',
        'results/miller-comparison/plots/distance-histogram/{sample}-{af}-subset1000.png',
        'results/miller-comparison/plots/distance-histogram/{sample}-{af}-full.png'
    params:
        out_prefix='results/miller-comparison/plots/distance-histogram/{sample}-{af}'
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=1,
        runtime_min=0,
        mem_total_mb=5*1024
    shell:
        'python scripts/plot-allele-distance-histogram.py -miller {input.miller} -vienna {input.vienna} -reference {input.reference} -af-cutoff {wildcards.af} -sample {wildcards.sample} -output {params.out_prefix}'

# histogram plot to compare the allele distance between VNTRs of mismatched sample vcfs
rule distance_histogram_cross:
    input:
        miller='results/miller-comparison/analysis/miller-tsv/{sample1}-sequence.tsv',
        vienna='results/miller-comparison/analysis/vienna-tsv/{sample2}-sequence.tsv',
        reference='results/miller-comparison/reference-vntrs.bed'
    output:
        'results/miller-comparison/plots/distance-histogram/{sample1}_{sample2}-{af}-subset100.png',
        'results/miller-comparison/plots/distance-histogram/{sample1}_{sample2}-{af}-subset1000.png',
        'results/miller-comparison/plots/distance-histogram/{sample1}_{sample2}-{af}-full.png'
    params:
        out_prefix='results/miller-comparison/plots/distance-histogram/{sample1}_{sample2}-{af}'
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=1,
        runtime_min=0,
        mem_total_mb=5*1024
    shell:
        'python scripts/plot-allele-distance-histogram.py -miller {input.miller} -vienna {input.vienna} -reference {input.reference} -af-cutoff {wildcards.af} -sample {wildcards.sample1}_{wildcards.sample2} -output {params.out_prefix}'

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

# creating VCF file for sample.
rule create_samplewise_vntr_vcf:
    input:
        sample_vntr='results/vamos-t2t/{sample}.vcf',
        ref_vntr='results/per-sample-bed/chm13-prep/chm13-vntr.vcf'
    output:
        'results/per-sample-bed/vcf/{sample}.vcf.gz'
    resources:
        runtime_hrs=1,
        runtime_min=0,
        mem_total_mb=1024
    shell:
        'python scripts/prepare-vntr-vcf.py -sample {input.sample_vntr} -ref {input.ref_vntr} | bgzip -c > {output}'

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