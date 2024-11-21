## Running QC using the HGSVC assemblies

# running VNTR calling on assemblies.
rule hgsvc_vntr_call:
    input:
        bam=config['path_to_hgsvc_alignments']+'{sample}.vrk-ps-sseq.asm-{haplotype}.t2tv2.sort.bam',
        sites='results/temp/vamos.T2T.processed.tsv'
    output:
        'results/hgsvc3-comparison/vntr-calls/{sample}-{haplotype}.vcf'
    log:
        'results/hgsvc3-comparison/vntr-calls/{sample}-{haplotype}.log'
    threads: 8
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=12,
        runtime_min=20,
        mem_total_mb=50*1024
    shell:
        'vamos --contig -b {input.bam} -r {input.sites} -s {wildcards.sample} -o {output} -t {threads}'

## VNTR Comparison Rules

# creating list of mismatched samples from the hgsvc list
mismatched_hgsvc_sample_set1 = []
mismatched_hgsvc_sample_set2 = []
mismatched_hgsvc_sample = []
for i in hgsvc_samples:
    for j in hgsvc_samples:
        if i == j:
            continue
        mismatched_hgsvc_sample_set1.append(i)
        mismatched_hgsvc_sample_set2.append(j)
        mismatched_hgsvc_sample.append(i+'_'+j)

# running VNTR comparison between HGSVC assemblies and ONT reads on the same sample
rule hgsvc_comp_same:
    input:
        ont='results/vamos-t2t/{sample}.vcf',
        hgsvc1='results/hgsvc3-comparison/vntr-calls/{sample}-hap1.vcf',
        hgsvc2='results/hgsvc3-comparison/vntr-calls/{sample}-hap2.vcf'
    output:
        'results/hgsvc3-comparison/qc/same-sample/{sample}-RU.svg',
        'results/hgsvc3-comparison/qc/same-sample/subset100-{sample}-RU.svg',
        'results/hgsvc3-comparison/qc/same-sample/{sample}.stats'
    log:
        'results/hgsvc3-comparison/qc/same-sample/{sample}.log'
    params:
        outprefix='results/hgsvc3-comparison/qc/same-sample/'
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=2*1024
    shell:
        'python scripts/plot-hgsvc-qc.py -hgsvc {input.hgsvc1},{input.hgsvc2} -ont {input.ont} -sample {wildcards.sample} -output {params.outprefix} 2> {log}'


# running VNTR comparison between HGSVC assemblies and ONT reads on the different sample
rule hgsvc_comp_diff:
    input:
        ont='results/vamos-t2t/{sample2}.vcf',
        hgsvc1='results/hgsvc3-comparison/vntr-calls/{sample1}-hap1.vcf',
        hgsvc2='results/hgsvc3-comparison/vntr-calls/{sample1}-hap2.vcf'
    output:
        'results/hgsvc3-comparison/qc/diff-sample/{sample1}_{sample2}-RU.svg',
        'results/hgsvc3-comparison/qc/diff-sample/subset100-{sample1}_{sample2}-RU.svg',
        'results/hgsvc3-comparison/qc/diff-sample/{sample1}_{sample2}.stats'
    log:
        'results/hgsvc3-comparison/qc/diff-sample/{sample1}_{sample2}.log'
    params:
        outprefix='results/hgsvc3-comparison/qc/diff-sample/'
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=2*1024
    shell:
        'python scripts/plot-hgsvc-qc.py -hgsvc {input.hgsvc1},{input.hgsvc2} -ont {input.ont} -sample {wildcards.sample1}_{wildcards.sample2} -output {params.outprefix} 2> {log}'

# creating the boxplot of PCCs
rule hgsvc_ru_pcc_boxplot:
    input:
        same_sample=expand('results/hgsvc3-comparison/qc/same-sample/{sample}.stats', sample=hgsvc_samples),
        diff_sample=expand('results/hgsvc3-comparison/qc/diff-sample/{diff_sample}.stats', diff_sample=mismatched_hgsvc_sample)
    output:
        'results/hgsvc3-comparison/qc/ru-boxplot-pcc.svg'
    params:
        same_sample=','.join(list(expand('results/hgsvc3-comparison/qc/same-sample/{sample}.stats', sample=hgsvc_samples))),
        diff_sample=','.join(list(expand('results/hgsvc3-comparison/qc/diff-sample/{diff_sample}.stats', diff_sample=mismatched_hgsvc_sample)))
    conda:
        '../envs/vamos.yml'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=2*1024
    shell:
        'python scripts/plot-ru-pcc-boxplot.py -same-sample {params.same_sample} -diff-sample {params.diff_sample} -output {output}'

# running stats for vntr-summary table on all the HGSVC assemblies (except HG00514. Warned not to use it by P. Ebert)
rule hgsvc_get_vntr_stats:
    input:
        hap1='results/hgsvc3-comparison/vntr-calls/{sample}-hap1.vcf',
        hap2='results/hgsvc3-comparison/vntr-calls/{sample}-hap2.vcf'
    output:
        'results/hgsvc3-comparison/vntr-calls/{sample}.stats'
    shell:
        'python scripts/get-vntr-stats-assemblies.py -hap1 {input.hap1} -hap2 {input.hap2} | awk \'$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n "}}\' > {output}'

# create the vntr summary table
rule hgsvc_vntr_summary_table:
    input:
        vcf=expand('results/hgsvc3-comparison/vntr-calls/{sample}.stats', sample=hgsvc_sample_list_all),
        sites='results/temp/vamos.T2T.processed.tsv'
    output:
        'results/hgsvc3-comparison/hgsvc-vntr-summary.bed'
    params:
        ','.join(list(expand('results/hgsvc3-comparison/vntr-calls/{sample}.stats', sample=hgsvc_sample_list_all)))
    conda:
        '../envs/comparison.yml'
    resources:
        runtime_hrs=5,
        runtime_min=0,
        mem_total_mb=96000
    shell:
        'python scripts/create-vntr-table.py -stats {params} -sites {input.sites} > {output}'

# plotting the vntr spread between hgsvc vntr calls (all assemblies) and ont vntr calls
rule plot_vntr_spread:
    input:
        hgsvc='results/hgsvc3-comparison/hgsvc-vntr-summary.bed',
        ont='results/vamos-t2t-summary.bed'
    output:
        'results/hgsvc3-comparison/vntr-spread/max-full.svg',
        'results/hgsvc3-comparison/vntr-spread/99-full.svg',
        'results/hgsvc3-comparison/vntr-spread/95-full.svg',
        'results/hgsvc3-comparison/vntr-spread/75-full.svg',
        'results/hgsvc3-comparison/vntr-spread/max-subset100.svg',
        'results/hgsvc3-comparison/vntr-spread/99-subset100.svg',
        'results/hgsvc3-comparison/vntr-spread/95-subset100.svg',
        'results/hgsvc3-comparison/vntr-spread/75-subset100.svg'
    params:
        par1='results/hgsvc3-comparison/vntr-spread/max',
        par2='results/hgsvc3-comparison/vntr-spread/99',
        par3='results/hgsvc3-comparison/vntr-spread/95',
        par4='results/hgsvc3-comparison/vntr-spread/75'
    conda:
        '../envs/vamos.yml'
    shell:
        '''
        python scripts/plot-vntr-spread.py -hgsvc {input.hgsvc} -ont {input.ont} -spread max -output {params.par1}
        python scripts/plot-vntr-spread.py -hgsvc {input.hgsvc} -ont {input.ont} -spread 99 -output {params.par2}
        python scripts/plot-vntr-spread.py -hgsvc {input.hgsvc} -ont {input.ont} -spread 95 -output {params.par3}
        python scripts/plot-vntr-spread.py -hgsvc {input.hgsvc} -ont {input.ont} -spread 75 -output {params.par4}
        '''