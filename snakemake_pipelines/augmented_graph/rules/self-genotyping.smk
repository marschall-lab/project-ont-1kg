### Self-genotyping Concordance with HG01258 ###

include: './get-sample-list.smk'

if config['pilot']:
    samples.sort()
    samples=samples[0:2]

wildcard_constraints:
    sample='|'.join(samples),
    regions='|'.join(['biallelic', 'multiallelic']),
    vartype='|'.join(['all', 'sv', 'indels', 'large-deletion', 'large-insertion', 'large-complex'])

# extract ground truth genotypes for sample
rule prepare_truth:
    input:
        "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/panel/giggles-ready_biallelic.vcf.gz"
    output:
        temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/truth/HG01258-truth-biallelic.vcf")
    conda:
        "../envs/basic.yml"
    priority: 1
    resources:
        mem_total_mb=20000
    log:
        "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/truth/HG01258-truth-biallelic.log"
    shell:
        "bcftools view --samples {wildcards.sample} {input} 2> {log} 1> {output}"

# plot alleles and create bed file with complex bubbles
rule alleles_per_bubble:
    input:
        "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/panel/giggles-ready_biallelic.vcf.gz"
    output:
        plot = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/alleles-per-bubble.pdf",
        bed = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/complex-bubbles.bed"
    conda:
        "../envs/basic.yml"
    resources:
        mem_total_mb=2000,
        runtime_hrs=1
    shell:
        "zcat {input} | python scripts/variant-statistics.py {output.plot} 1 > {output.bed}"

# prepare beds for biallelic and complex graph regions
rule prepare_beds:
    input:
        bed = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/complex-bubbles.bed",
        fai = '/gpfs/project/projects/medbioinf/users/spani/files/fasta/HPRC/reference/chm13v2.0_maskedY_rCRS.fa.gz.fai'
    output:
        bed = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/biallelic-bubbles.bed",
        tmp = temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/biallelic-bubbles.fai")
    conda:
        "../envs/basic.yml"
    shell:
        """
        sort -k1,1d -k 2,2n -k 3,3n {input.fai} > {output.tmp}
        bedtools complement -i {input.bed} -g {output.tmp} > {output.bed}
        """

# create empty list of untypable variant IDs
rule create_empty_list_of_untypable_ids:
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/HG01258-untypable-ids.tsv')
    shell:
        'touch {output}'

# determine untypable IDs
rule remove_untypable_callset:
    input:
        vcf = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/genotypes/HG01258-biallelic.vcf.gz",
        tmp = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/HG01258-untypable-ids.tsv"
    output:
        vcf = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/callset/HG01258-typable-{vartype}.vcf.gz",
        tbi = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/callset/HG01258-typable-{vartype}.vcf.gz.tbi"
    resources:
        mem_total_mb = 200,
        runtime_hrs = 1
    conda:
        "../envs/basic.yml"
    shell:
        """
        zcat {input.vcf} | python scripts/skip-untypable.py {input.tmp} | python scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# determine untypable IDs
rule remove_untypable_truth:
    input:
        vcf = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/truth/HG01258-truth-biallelic.vcf.gz",
        tmp = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/HG01258-untypable-ids.tsv"
    output:
        vcf = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/truth/HG01258-truth-biallelic-{vartype}.vcf.gz",
        tbi = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/truth/HG01258-truth-biallelic-{vartype}.vcf.gz.tbi"
    resources:
        mem_total_mb = 200,
        runtime_hrs = 1
    conda:
        "../envs/basic.yml"
    shell:
        """
        zcat {input.vcf} | python scripts/skip-untypable.py {input.tmp} | python scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

def region_to_bed(wildcards):
    if wildcards.regions == "biallelic":
        return "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/biallelic-bubbles.bed".format(callset=wildcards.callset)
    if wildcards.regions == "multiallelic":
        return "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/complex-bubbles.bed".format(callset=wildcards.callset)
    assert(False)

# determine the variants that went into re-typing per category
rule collect_typed_variants:
    input:
        callset = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/panel/giggles-ready_biallelic.vcf.gz",
        regions= region_to_bed,
        ids="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/HG01258-untypable-ids.tsv"
    output:
        "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/genotyped-ids/HG01258-{regions}-{vartype}.tsv"
    conda:
        "../envs/basic.yml"
    resources:
        mem_total_mb=5000
    shell:
        "zcat {input.callset} | python scripts/skip-untypable.py {input.ids} | python scripts/extract-varianttype.py {wildcards.vartype} | bedtools intersect -header -a - -b {input.regions} -u -f 0.5 | python scripts/get_ids.py > {output}"

# compute concordances
rule genotype_concordances:
    input:
        callset = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/callset/HG01258-typable-{vartype}.vcf.gz",
        baseline = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/truth/HG01258-truth-biallelic-{vartype}.vcf.gz",
        regions = region_to_bed,
        typed_ids = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/genotyped-ids/HG01258-{regions}-{vartype}.tsv"
    output:
        tmp_vcf1 = temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/concordance/{regions}_{vartype}_base.vcf"),
        tmp_vcf2 = temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/concordance/{regions}_{vartype}_call.vcf"),
        summary = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/concordance/{regions}_{vartype}/summary.txt"
    conda:
        "../envs/basic.yml"
    log:
        "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/concordance/{regions}_{vartype}/summary.log"
    resources:
        mem_total_mb = 400,
        runtime_hrs = 0,
        runtime_min = 40
    shell:
        """
        bedtools intersect -header -a {input.baseline} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf1}
        bedtools intersect -header -a {input.callset} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf2}
        python scripts/genotype-evaluation.py {output.tmp_vcf1} {output.tmp_vcf2} {input.typed_ids} --qual 0 2> {log} 1> {output.summary}
        """

# summarize concordance results
rule summarize_concordance:
    input:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/self-genotyping/concordance/{regions}_{vartype}/summary.txt', regions=['biallelic', 'multiallelic'],vartype=['all', 'sv', 'indels', 'large-deletion', 'large-insertion', 'large-complex'])
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/concordance/summary.tsv'
    params:
        dir='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/concordance/',
        vartype=','.join(['all', 'sv', 'indels', 'large-deletion', 'large-insertion', 'large-complex']),
        regions='biallelic,multiallelic'
    conda:
        "../envs/basic.yml"
    shell:
        'python scripts/concordance-summary.py {params.dir} {params.regions} {params.vartype} > {output}'