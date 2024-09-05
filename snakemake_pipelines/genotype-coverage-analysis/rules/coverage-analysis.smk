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
        vcf=config['path_to_phased_vcf'],
        tsv='results/coverage-experiments/sample-by-coverage_{cov_range}.tsv'
    output:
        vcf=temp('results/coverage-experiments/subsampled-vcf/{cov_range}.vcf')
    conda:
        '../envs/coverage-analysis.yml'
    shell:
        'bcftools view --force-samples -S {input.tsv} {input.vcf} > {output.vcf}'

# creating the QC table
rule create_qc_table:
    input:
        vcf='results/coverage-experiments/subsampled-vcf/{cov_range}.vcf',
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
        'python scripts/qc-table.py -bi-panel {input.panel_bi} -bi-callset {input.vcf} -multi-panel {input.panel_multi} -sample-sheet {input.sample_sheet}'

# creating a checkpoint with all the TSV files that have been generated
checkpoint aggregate_tsvs:
    input:
        expand('results/coverage-experiments/qc-tables/{cov_range}.chk', cov_range=coverage_ranges)
    output:
        directory('results/coverage-experiments/qc-tables-selected/')
    shell:
        'cp results/coverage-experiments/qc-tables/*.tsv results/coverage-experiments/qc-tables-selected'

# creating the plots of the selected tsvs
rule make_plots:
    input:
        'results/coverage-experiments/qc-tables-selected/{cov_ranges}.tsv'
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
    
# aggregating all the pdfs generated (only for the all category of sv types)
def aggregate_pdfs(wildcards):
    chk_out=checkpoints.aggregate_tsvs.get(**wildcards).output[0]
    outfiles=[i for i in expand('results/coverage-experiments/plots/{cov_range}/all.pdf', cov_range=glob_wildcards(os.path.join(chk_out, "{cov_range}.tsv")).cov_range)]
    return outfiles

# making a single pdf with all the images
rule plot_to_one_pdf:
    input:
        aggregate_pdfs
    output:
        'results/coverage-experiments/plots/all.pdf'
    shell:
        'pdfunite {input} {output}'