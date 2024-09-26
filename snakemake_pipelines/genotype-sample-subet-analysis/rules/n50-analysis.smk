# This pipeline is for separating samples based on n50 and analysing the n50 ranges.

# list of parameters with their values
binsize = 10
max_n50 = 60      # hardcoded maximum n50 to create n50 ranges
n50_ranges = [str(n50)+'-'+str(n50+binsize-1) for n50 in range(0, max_n50, binsize)]
sv_types = ['large_deletions', 'large_insertions', 'large_complex', 'all']

# global wildcard constraints
wildcard_constraints:
    n50_range='|'.join(n50_ranges),
    sv_type='|'.join(sv_types)

# process sample data and bin by n50
rule bin_sample_data_n50:
    input:
        'resources/1k_ont_data_overview-samples_and_sequencing.tsv'
    output:
        expand('results/n50-experiments/sample-by-n50_{n50_range}.tsv', n50_range=n50_ranges)
    params:
        outprefix='results/n50-experiments/sample-by-n50',
        bs=binsize
    conda:
        '../envs/analysis.yml'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=5*1024
    shell:
        'python scripts/bin-sample-data-by-n50.py -binsize {params.bs} -tsv {input} -output {params.outprefix}'

# subsetting vcf to samples based on n50
rule subset_vcf_n50:
    input:
        callset_vcf=config['path_to_callset_vcf'],
        sniffles_vcf=config['path_to_sniffles_vcf'],
        delly_vcf=config['path_to_delly_vcf'],
        svarp_vcf=config['path_to_svarp_vcf'],
        tsv='results/n50-experiments/sample-by-n50_{n50_range}.tsv'
    output:
        temp('results/n50-experiments/subsampled-vcf/{n50_range}-callset.vcf'),
        temp('results/n50-experiments/subsampled-vcf/{n50_range}-sniffles.vcf'),
        temp('results/n50-experiments/subsampled-vcf/{n50_range}-delly.vcf'),
        temp('results/n50-experiments/subsampled-vcf/{n50_range}-svarp.vcf')
    params:
        outdir='results/n50-experiments/subsampled-vcf/{n50_range}'
    conda:
        '../envs/analysis.yml'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=3*1024
    shell:
        '''
        bcftools view --force-samples -S {input.tsv} {input.callset_vcf} > {params.outdir}-callset.vcf
        bcftools view --force-samples -S {input.tsv} {input.sniffles_vcf} > {params.outdir}-sniffles.vcf
        bcftools view --force-samples -S {input.tsv} {input.delly_vcf} > {params.outdir}-delly.vcf
        bcftools view --force-samples -S {input.tsv} {input.svarp_vcf} > {params.outdir}-svarp.vcf
        '''

# creating the QC table
rule create_qc_table_n50:
    input:
        vcf='results/n50-experiments/subsampled-vcf/{n50_range}-callset.vcf',
        panel_bi=config['path_to_panel_bi'],
        panel_multi=config['path_to_panel_multi'],
        sample_sheet='resources/sample_sheet.tsv'
    output:
        'results/n50-experiments/qc-tables/{n50_range}.chk'
    params:
        outprefix='results/n50-experiments/qc-tables/{n50_range}'
    log:
        'results/n50-experiments/qc-tables/{n50_range}.log'
    conda:
        '../envs/analysis.yml'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=5*1024
    shell:
        'python scripts/qc-table.py -bi-panel {input.panel_bi} -bi-callset {input.vcf} -multi-panel {input.panel_multi} -sample-sheet {input.sample_sheet} -output {params.outprefix}'

# creating a checkpoint with all the TSV files that have been generated
checkpoint aggregate_tsvs_n50:
    input:
        expand('results/n50-experiments/qc-tables/{n50_range}.chk', n50_range=n50_ranges)
    output:
        directory('results/n50-experiments/qc-tables-selected/')
    shell:
        'mkdir -p {output} && cp results/n50-experiments/qc-tables/*.tsv {output}'

# creating the plots of the selected tsvs
rule make_plots_n50:
    input:
        'results/n50-experiments/qc-tables-selected/{n50_range}.tsv'
    output:
        expand('results/n50-experiments/qc-plots/{{n50_range}}/{sv_type}.pdf', sv_type=sv_types)
    params:
        out='results/n50-experiments/qc-plots/{n50_range}/'
    conda:
        '../envs/analysis.yml'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=5*1024
    shell:
        'python scripts/plot-hwe-af.py -table {input} -output {params.out}'
    
# aggregating all the relevant n50 ranges
def get_ranges_n50(wildcards):
    chk_out=checkpoints.aggregate_tsvs_n50.get(**wildcards).output[0]
    RNG = glob_wildcards(os.path.join(chk_out, "{rng}.tsv"))
    return RNG.rng

def aggregate_pdfs_n50(wildcards):
    ranges=get_ranges_n50(wildcards)
    return expand('results/n50-experiments/qc-plots/{n50_range}/all.pdf', n50_range=ranges)

# making a single pdf with all the images
rule plot_to_one_pdf_n50:
    input:
        aggregate_pdfs_n50
    output:
        'results/n50-experiments/qc-plots/all.pdf'
    conda:
        '../envs/analysis.yml'
    shell:
        'pdfcombine --force -o {output} {input}'

# process the sub-sampled vcfs to get the number of SVs per sample
rule sv_number_per_sample_n50:
    input:
        callset='results/n50-experiments/subsampled-vcf/{n50_range}-callset.vcf',
        sniffles='results/n50-experiments/subsampled-vcf/{n50_range}-sniffles.vcf',
        delly='results/n50-experiments/subsampled-vcf/{n50_range}-delly.vcf',
        svarp='results/n50-experiments/subsampled-vcf/{n50_range}-svarp.vcf'
    output:
        'results/n50-experiments/sv_count_per_sample/{n50_range}.tsv'
    conda:
        '../envs/analysis.yml'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=5*1024
    shell:
        'python scripts/count-svs-per-sample.py -callset {input.callset} -sniffles {input.sniffles} -delly {input.delly} -svarp {input.svarp} > {output}'


def aggregate_sv_count_tsvs_n50(wildcards):
    ranges=get_ranges_n50(wildcards)
    return expand('results/n50-experiments/sv_count_per_sample/{n50_range}.tsv', n50_range=ranges)

# plot all the sv numbers for various n50 ranges
rule plot_sv_counts_n50:
    input:
        aggregate_sv_count_tsvs_n50
    output:
        'results/n50-experiments/sv_count_per_sample/plot.svg'
    params:
        lambda wildcards: ','.join([i for i in aggregate_sv_count_tsvs_n50(wildcards)])
    conda:
        '../envs/analysis.yml'
    shell:
        'python scripts/plot-sv-counts.py -tsvs {params} -output {output}'
