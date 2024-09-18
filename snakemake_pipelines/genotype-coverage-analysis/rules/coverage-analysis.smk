# This pipeline is for separating samples based on coverage and analysing the coverage ranges.

# list of parameters with their values
binsize = 10
max_coverage = 140      # hardcoded maximum coverage to create coverage ranges
coverage_ranges = [str(cov)+'-'+str(cov+binsize-1) for cov in range(0, max_coverage, binsize)]
sv_types = ['large_deletions', 'large_insertions', 'large_complex', 'all']

# global wildcard constraints
wildcard_constraints:
    cov_ranges='|'.join(coverage_ranges),
    sv_type='|'.join(sv_types)

# process sample data and bin by coverage
rule bin_sample_data:
    input:
        'resources/1k_ont_data_overview-samples_and_sequencing.tsv'
    output:
        expand('results/coverage-experiments/sample-by-coverage_{cov_range}.tsv', cov_range=coverage_ranges)
    params:
        outprefix='results/coverage-experiments/sample-by-coverage',
        bs=binsize
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=5*1024
    shell:
        'python scripts/bin-sample-data-by-coverage.py -binsize {params.bs} -tsv {input} -output {params.outprefix}'

# subsetting vcf to samples based on coverage
rule subset_vcf:
    input:
        callset_vcf=config['path_to_phased_vcf'],
        sniffles_vcf=config['path_to_sniffles_vcf'],
        delly_vcf=config['path_to_delly_vcf'],
        svarp_vcf=config['path_to_svarp_vcf'],
        tsv='results/coverage-experiments/sample-by-coverage_{cov_range}.tsv'
    output:
        temp('results/coverage-experiments/subsampled-vcf/{cov_range}-callset.vcf'),
        temp('results/coverage-experiments/subsampled-vcf/{cov_range}-sniffles.vcf'),
        temp('results/coverage-experiments/subsampled-vcf/{cov_range}-delly.vcf'),
        temp('results/coverage-experiments/subsampled-vcf/{cov_range}-svarp.vcf')
    params:
        outdir='results/coverage-experiments/subsampled-vcf/{cov_range}'
    conda:
        '../envs/coverage-analysis.yml'
    shell:
        '''
        bcftools view --force-samples -S {input.tsv} {input.callset_vcf} > {params.outdir}-callset.vcf
        bcftools view --force-samples -S {input.tsv} {input.sniffles_vcf} > {params.outdir}-sniffles.vcf
        bcftools view --force-samples -S {input.tsv} {input.delly_vcf} > {params.outdir}-delly.vcf
        bcftools view --force-samples -S {input.tsv} {input.svarp_vcf} > {params.outdir}-svarp.vcf
        '''

# creating the QC table
rule create_qc_table:
    input:
        vcf='results/coverage-experiments/subsampled-vcf/{cov_range}-callset.vcf',
        panel_bi=config['path_to_panel_bi'],
        panel_multi=config['path_to_panel_multi'],
        sample_sheet='resources/sample_sheet.tsv'
    output:
        'results/coverage-experiments/qc-tables/{cov_range}.chk'
    params:
        outprefix='results/coverage-experiments/qc-tables/{cov_range}'
    log:
        'results/coverage-experiments/qc-tables/{cov_range}.log'
    conda:
        '../envs/coverage-analysis.yml'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=5*1024
    shell:
        'python scripts/qc-table.py -bi-panel {input.panel_bi} -bi-callset {input.vcf} -multi-panel {input.panel_multi} -sample-sheet {input.sample_sheet} -output {params.outprefix}'

# creating a checkpoint with all the TSV files that have been generated
checkpoint aggregate_tsvs:
    input:
        expand('results/coverage-experiments/qc-tables/{cov_range}.chk', cov_range=coverage_ranges)
    output:
        directory('results/coverage-experiments/qc-tables-selected/')
    shell:
        'mkdir -p {output} && cp results/coverage-experiments/qc-tables/*.tsv {output}'

# creating the plots of the selected tsvs
rule make_plots:
    input:
        'results/coverage-experiments/qc-tables-selected/{cov_range}.tsv'
    output:
        expand('results/coverage-experiments/plots/{{cov_range}}/{sv_type}.pdf', sv_type=sv_types)
    params:
        out='results/coverage-experiments/plots/{cov_range}/'
    conda:
        '../envs/coverage-analysis.yml'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=5*1024
    shell:
        'python scripts/plot-hwe-af.py -table {input} -output {params.out}'
    
# aggregating all the relevant coverage ranges
def get_ranges(wildcards):
    chk_out=checkpoints.aggregate_tsvs.get(**wildcards).output[0]
    RNG = glob_wildcards(os.path.join(chk_out, "{rng}.tsv"))
    return RNG.rng

def aggregate_pdfs(wildcards):
    ranges=get_ranges(wildcards)
    return expand('results/coverage-experiments/plots/{cov_range}/all.pdf', cov_range=ranges)

# making a single pdf with all the images
rule plot_to_one_pdf:
    input:
        aggregate_pdfs
    output:
        'results/coverage-experiments/plots/all.pdf'
    shell:
        'pdfcombine --force -o {output} {input}'

# process the sub-sampled vcfs to get the number of SVs per sample
rule sv_number_per_sample:
    input:
        callset='results/coverage-experiments/subsampled-vcf/{cov_range}-callset.vcf',
        sniffles='results/coverage-experiments/subsampled-vcf/{cov_range}-sniffles.vcf',
        delly='results/coverage-experiments/subsampled-vcf/{cov_range}-delly.vcf',
        svarp='results/coverage-experiments/subsampled-vcf/{cov_range}-svarp.vcf'
    output:
        'results/sv_count_per_sample/{cov_range}.tsv'
    conda:
        '../envs/coverage-analysis.yml'
    shell:
        'python scripts/count-svs-per-sample.py -callset {input.callset} -sniffles {input.sniffles} -delly {input.delly} -svarp {input.svarp} -range {wildcards.range} > {output}'


def aggregate_sv_count_tsvs(wildcards):
    ranges=get_ranges(wildcards)
    return expand('results/sv_count_per_sample/{cov_range}.tsv', cov_range=ranges)

# plot all the sv numbers for various coverage ranges
rule plot_sv_counts:
    input:
        aggregate_sv_count_tsvs
    output:
        'results/sv_count_per_sample/plot.pdf'
    params:
        lambda wildcards: ','.join([i for i in aggregate_sv_count_tsvs(wildcards)])
    conda:
        '../envs/coverage-analysis.yml'
    shell:
        'python scripts/plot-sv-counts.py -tsvs {params} -output {output}'
