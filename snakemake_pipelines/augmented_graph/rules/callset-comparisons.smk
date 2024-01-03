# callsets = [c for c in config["callset_vcfs"].keys()]
max_af = 1
min_af = 0

# basic counting stats of the vcfs
rule vcf_stats:
    input:
        lambda wildcards: config["callset_vcfs"][wildcards.callset]
    output:
        "results/vcf-stats/{callset}.stats"
    conda:
        "../envs/comparison.yml"
    wildcard_constraints:
        callset = "|".join(callsets)
    shell:
        "zcat {input} | python scripts/set-pass.py | bcftools view -f PASS --min-af {min_af} --max-af {max_af} | python3 workflow/scripts/vcf_stats.py > {output}"

# add tag of variant types
rule add_tags:
    input:
        lambda wildcards: config["callset_vcfs"][wildcards.callset]
    output:
        "results/vcfs/{callset}-tagged.vcf"
    wildcard_constraints:
        callset = "|".join(callsets)
    conda:
        "../envs/comparison.yml"
    shell:
        "zcat {input} | python scripts/set-pass.py | bcftools view -f PASS --min-af {min_af} --max-af {max_af} | python3 workflow/scripts/add-svtags.py > {output}"

# extract sample which will be compared
rule extract_sample:
    input:
        "results/vcfs/{callset}-tagged.vcf"
    output:
        "results/vcfs/{callset}-tagged-{sample}.vcf"
    conda:
        "../envs/comparison.yml"
    shell:
        "bcftools view --samples {wildcards.sample} {input} | bcftools view --min-ac 1 > {output}"


def intersect_vcfs_files(wildcards):
    files = []
    for c in callsets:
        if wildcards.sample == "all":
            files.append("results/vcfs/{callset}-tagged.vcf".format(callset = c))
        else:
            files.append("results/vcfs/{callset}-tagged-{sample}.vcf".format(callset = c, sample = wildcards.sample))
    return files

# intersecting the vcf files and creating upset plots
rule intersect_vcfs:
    input:
        intersect_vcfs_files
    output:
        tsv="results/intersection/intersection_{sample}.tsv",
        vcf="results/intersection/intersection_{sample}.vcf",
        pdf="results/intersection/intersection_{sample}.pdf",
        plot="results/intersection/intersection_upset_{sample}.pdf"
    conda:
        "../envs/plotting.yml"
    log:
        intersect="results/intersection/intersection_{sample}.log",
        plot="results/intersection/plotting_{sample}.log"
    params:
        names = callsets,
        columns = ["in_" + c for c in callsets]
    shell:
        """
        python scripts/intersect_callsets.py intersect -c {input} -n {params.names} -t {output.tsv} -v {output.vcf} -p {output.pdf} &> {log.intersect}
        python scripts/plot-comparison-upset.py -t {output.tsv} -o {output.plot} -n {params.columns} &> {log.plot}
        """