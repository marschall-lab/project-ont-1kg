# sort miller vcf
rule sort_miller_vcf:
    input:
        config['path_to_miller_vcf']
    output:
        temp('results/analysis/vamos-miller.sorted.vcf')
    shell:
        "zcat {input} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}"

# create TSV of miller vcf
rule parse_miller_vcf:
    input:
        'results/analysis/vamos-miller.sorted.vcf'
    output:
        'results/analysis/miller-tsv/{sample}-count.tsv',
        'results/analysis/miller-tsv/{sample}-sequence.tsv'
    params:
        out_prefix='results/analysis/miller-tsv/{sample}'
    shell:
        'python scripts/vamos-parse-miller.py -vcf {input} -sample {wildcards.sample} -output {params.out_prefix}'

# create TSV of vienna vcf
rule parse_vienna_vcf:
    input:
        'results/vamos-results-vienna/{sample}.vcf'
    output:
        'results/analysis/vienna-tsv/{sample}-count.tsv',
        'results/analysis/vienna-tsv/{sample}-sequence.tsv'
    params:
        out_prefix='results/analysis/vienna-tsv/{sample}'
    shell:
        'python scripts/vamos-parse-vienna.py -vcf {input} -sample {wildcards.sample} -output {params.out_prefix}'

# density plot to compare the count of repeating units in both vcfs
rule ru_densityplot:
    input:
        miller='results/analysis/miller-tsv/{sample}-count.tsv',
        vienna='results/analysis/vienna-tsv/{sample}-count.tsv'
    output:
        'results/plots/ru-density/{sample}.png'
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
        miller='results/analysis/miller-tsv/{sample1}-count.tsv',
        vienna='results/analysis/vienna-tsv/{sample2}-count.tsv'
    output:
        'results/plots/ru-density/{sample1}_{sample2}.png'
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
        miller='results/analysis/miller-tsv/{sample}-sequence.tsv',
        vienna='results/analysis/vienna-tsv/{sample}-sequence.tsv',
        reference='results/reference-vntrs.bed'
    output:
        'results/plots/distance-histogram/{sample}-{af}-subset1000.png',
        'results/plots/distance-histogram/{sample}-{af}-full.png'
    params:
        out_prefix='results/plots/distance-histogram/{sample}-{af}'
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
        miller='results/analysis/miller-tsv/{sample1}-sequence.tsv',
        vienna='results/analysis/vienna-tsv/{sample2}-sequence.tsv',
        reference='results/reference-vntrs.bed'
    output:
        'results/plots/distance-histogram/{sample1}_{sample2}-{af}-subset1000.png',
        'results/plots/distance-histogram/{sample1}_{sample2}-{af}-full.png'
    params:
        out_prefix='results/plots/distance-histogram/{sample1}_{sample2}-{af}'
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=1,
        runtime_min=0,
        mem_total_mb=5*1024
    shell:
        'python scripts/plot-allele-distance-histogram.py -miller {input.miller} -vienna {input.vienna} -reference {input.reference} -af-cutoff {wildcards.af} -sample {wildcards.sample1}_{wildcards.sample2} -output {params.out_prefix}'