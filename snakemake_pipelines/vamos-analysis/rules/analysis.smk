# sort miller vcf
rule sort_miller_vcf:
    input:
        config['path_to_miller_vcf']
    output:
        temp('results/analysis/vamos-miller.sorted.vcf')
    shell:
        "zcat {input} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}"

# create TSV of repeating unit counts of miller vcf
rule ru_count_miller:
    input:
        'results/analysis/vamos-miller.sorted.vcf'
    output:
        'reulsts/analysis/miller-tsv/{sample}.tsv'
    shell:
        'python scripts/vamos-parse-miller.py -vcf {input} -sample {wildcards.sample} > {output}'

# create TSV of repeating unit counts of vienna vcf
rule ru_count_vienna:
    input:
        'results/vamos-results-vienna/{sample}.vcf'
    output:
        'results/analysis/vienna-tsv/{sample}.tsv'
    shell:
        'python scripts/vamos-parse-vienna.py -vcf {input} -sample {wildcards.sample} > {output}'

# scatter plot to compare the count of repeating units in both vcfs
rule ru_scatterplot:
    input:
        miller='reulsts/analysis/miller-tsv/{sample}.tsv',
        vienna='results/analysis/vienna-tsv/{sample}.tsv'
    output:
        'results/plots/ru-scatter/{sample}.png'
    conda:
        '../envs/vamos.yml'
    shell:
        'python scripts/plot-ru-count-scatter.py -miller {input.miller} -vienna {input.vienna} -output {output}'