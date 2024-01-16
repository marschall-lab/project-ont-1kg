### Self-genotyping Concordance with HG01258 ###

include: './get-sample-list.smk'

if config['pilot']:
    samples.sort()
    samples=samples[0:2]

wildcard_constraints:
    sample='|'.join(samples),
    regions='|'.join(['biallelic', 'multiallelic']),
    vartype='|'.join(['all', 'sv', 'large-deletion', 'large-insertion', 'large-complex', 'large-deletion-1000', 'large-insertion-1000', 'large-complex-1000', 'large-deletion-10000', 'large-insertion-10000', 'large-complex-10000', 'large-deletion-50000', 'large-insertion-50000', 'large-complex-50000', 'large-deletion-above50000', 'large-insertion-above50000', 'large-complex-above50000']),
    representation='biallelic|multiallelic'


# subset vcfs based on allele frequency
rule subset_vcf_af:
    input:
        panel='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/panel/giggles-ready_biallelic.vcf.gz',
        callset='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/genotypes/HG01258-biallelic.vcf.gz'
    output:
        panel='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/truth/subset/giggles-ready_biallelic_{min_af}-{max_af}.vcf',
        callset='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/callset/subset/HG01258-biallelic_{min_af}-{max_af}.vcf'
    conda:
        '../envs/basic.yml'
    shell:
        '''
        python scripts/subset-by-af.py -panel {input.panel} -callset {input.callset} -mode panel -min-af {wildcards.min_af} -max-af {wildcards.max_af} -output-panel {output.panel} -output-callset {output.callset}
        '''

# extract ground truth genotypes for sample
rule prepare_truth:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/truth/subset/giggles-ready_biallelic_{min_af}-{max_af}.vcf'
    output:
        temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/truth/HG01258-truth-biallelic_{min_af}-{max_af}.vcf")
    conda:
        "../envs/basic.yml"
    resources:
        mem_total_mb=20000
    shell:
        "bcftools view --samples HG01258 {input} > {output}"

# plot alleles and create bed file with complex bubbles
rule alleles_per_bubble:
    input:
        "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/panel/giggles-ready_{representation}.vcf.gz"
    output:
        plot = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/{representation}-alleles-per-bubble.pdf",
        bed = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/{representation}-complex-bubbles.bed"
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
        bed = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/{representation}-complex-bubbles.bed",
        fai = '/gpfs/project/projects/medbioinf/users/spani/files/fasta/HPRC/reference/chm13v2.0_maskedY_rCRS.fa.gz.fai'
    output:
        bed = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/{representation}-biallelic-bubbles.bed",
        tmp = temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/{representation}-biallelic-bubbles.fai")
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
        vcf = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/callset/subset/HG01258-biallelic_{min_af}-{max_af}.vcf",
        tmp = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/HG01258-untypable-ids.tsv"
    output:
        vcf = temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/callset/HG01258-callset-typable_{min_af}-{max_af}_{vartype}.vcf.gz"),
        tbi = temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/callset/HG01258-callset-typable_{min_af}-{max_af}_{vartype}.vcf.gz.tbi")
    resources:
        mem_total_mb = 1000,
        runtime_hrs = 1
    conda:
        "../envs/basic.yml"
    shell:
        """
        cat {input.vcf} | python scripts/skip-untypable.py {input.tmp} | python scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# determine untypable IDs
rule remove_untypable_truth:
    input:
        vcf = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/truth/HG01258-truth-biallelic_{min_af}-{max_af}.vcf",
        tmp = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/HG01258-untypable-ids.tsv"
    output:
        vcf = temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/truth/HG01258-truth-typable_{min_af}-{max_af}_{vartype}.vcf.gz"),
        tbi = temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/truth/HG01258-truth-typable_{min_af}-{max_af}_{vartype}.vcf.gz.tbi")
    resources:
        mem_total_mb = 1000,
        runtime_hrs = 1
    conda:
        "../envs/basic.yml"
    shell:
        """
        cat {input.vcf} | python scripts/skip-untypable.py {input.tmp} | python scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

def region_to_bed(wildcards):
    if wildcards.regions == "biallelic":
        return "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/{representation}-biallelic-bubbles.bed".format(callset=wildcards.callset, representation=wildcards.representation)
    if wildcards.regions == "multiallelic":
        return "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/{representation}-complex-bubbles.bed".format(callset=wildcards.callset, representation=wildcards.representation)
    assert(False)

# determine the variants that went into re-typing per category
rule collect_typed_variants:
    input:
        callset = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/truth/subset/giggles-ready_biallelic_{min_af}-{max_af}.vcf",
        regions= region_to_bed,
        ids="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/HG01258-untypable-ids.tsv"
    output:
        "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/genotyped-ids-{representation}/HG01258-{regions}-{vartype}_{min_af}-{max_af}.tsv"
    conda:
        "../envs/basic.yml"
    resources:
        mem_total_mb=5000
    shell:
        "cat {input.callset} | python scripts/skip-untypable.py {input.ids} | python scripts/extract-varianttype.py {wildcards.vartype} | bedtools intersect -header -a - -b {input.regions} -u -f 0.5 | python scripts/get_ids.py > {output}"

# compute concordances
rule genotype_concordances:
    input:
        callset = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/callset/HG01258-callset-typable_{min_af}-{max_af}_{vartype}.vcf.gz",
        baseline = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/truth/HG01258-truth-typable_{min_af}-{max_af}_{vartype}.vcf.gz",
        regions = region_to_bed,
        typed_ids = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/genotyped-ids-{representation}/HG01258-{regions}-{vartype}_{min_af}-{max_af}.tsv"
    output:
        tmp_vcf1 = temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/concordance-{representation}/{regions}_{vartype}_{min_af}-{max_af}_base.vcf"),
        tmp_vcf2 = temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/concordance-{representation}/{regions}_{vartype}_{min_af}-{max_af}_call.vcf"),
        summary = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/concordance-{representation}/{regions}_{vartype}_{min_af}-{max_af}/summary.txt"
    conda:
        "../envs/basic.yml"
    log:
        "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/concordance-{representation}/{regions}_{vartype}_{min_af}-{max_af}/summary.log"
    resources:
        mem_total_mb = 2000,
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
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/self-genotyping/concordance-{{representation}}/{regions}_{vartype}_{{min_af}}-{{max_af}}/summary.txt', regions=['biallelic', 'multiallelic'],vartype=['all', 'sv', 'large-deletion', 'large-insertion', 'large-complex', 'large-deletion-1000', 'large-insertion-1000', 'large-complex-1000', 'large-deletion-10000', 'large-insertion-10000', 'large-complex-10000', 'large-deletion-50000', 'large-insertion-50000', 'large-complex-50000', 'large-deletion-above50000', 'large-insertion-above50000', 'large-complex-above50000'])
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/concordance-{representation}/summary_{min_af}-{max_af}.tsv'
    params:
        dir='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/self-genotyping/concordance-{representation}/',
        vartype=','.join(['all', 'sv', 'large-deletion', 'large-insertion', 'large-complex', 'large-deletion-1000', 'large-insertion-1000', 'large-complex-1000', 'large-deletion-10000', 'large-insertion-10000', 'large-complex-10000', 'large-deletion-50000', 'large-insertion-50000', 'large-complex-50000', 'large-deletion-above50000', 'large-insertion-above50000', 'large-complex-above50000']),
        regions='biallelic,multiallelic'
    conda:
        "../envs/basic.yml"
    shell:
        'python scripts/concordance-summary.py {params.dir} {params.regions} {params.vartype} {wildcards.min_af} {wildcards.max_af} > {output}'