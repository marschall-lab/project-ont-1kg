include: './get-sample-list.smk'

if config['pilot']:
    samples.sort()
    samples=samples[0:2]

wildcard_constraints:
    sample='|'.join(samples),
    pop='|'.join(['AFR', 'AMR', 'EAS', 'EUR', 'SAS']),
    sv_type='|'.join(['COMPLEX', 'DEL', 'INS']),
    ranges='|'.join(['0-1kbp', '1kbp-10kbp', '10kbp-100kbp', '100kbp-1Mbp']),
    regions='|'.join(['biallelic', 'multiallelic']),
    vartype='|'.join(['all', 'sv', 'indels', 'large-deletion', 'large-insertion', 'large-complex'])

# run genotyping
rule giggles:
    input:
        reads='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/data/fasta/{sample}.fasta.gz',
        fai='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/data/fasta/{sample}.fasta.gz.fai',
        gzi='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/data/fasta/{sample}.fasta.gz.gzi',
        haplotag='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/haplotagging-results/GRCh38/{sample}/{sample}.tsv',
        alignment='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/ont-alignments/{sample}.sorted.gaf.gz',
        alignment_index='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/ont-alignments/{sample}.sorted.gaf.gz.gai',
        gfa='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa',
        vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/panel/giggles-ready_multiallelic.vcf.gz'
    output:
        vcf=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/genotypes/{sample}-multiallelic.vcf')
    log:
        stderr='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/genotypes/{sample}-multiallelic.stderr'
    resources:
        mem_total_mb=20000,
        runtime_hrs=71,
        runtime_min=59
    shell:
        """
        set +u
        source ~/.bashrc
        conda activate giggles-dev
        set -u
        giggles --version > {log}
        giggles genotype --read-fasta {input.reads} --sample {wildcards.sample} --realign-mode edit -o {output.vcf} --rgfa {input.gfa} --haplotag-tsv {input.haplotag} {input.vcf} {input.alignment} 2>> {log.stderr}
        """


# get the biallelic and multiallelic VCFs for the genotypes
rule convert_sample_vcf_biallelic:
    input:
        sample_vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/genotypes/{sample}-multiallelic.vcf',
        biallelic_vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/panel/giggles-ready_biallelic.vcf.gz'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/genotypes/{sample}-biallelic.vcf')
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=3000
    shell:
        "cat {input.sample_vcf} | python scripts/convert-to-biallelic.py {input.biallelic_vcf} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}"


# merge multiallelic records to give multisample vcf
rule merge_vcf_to_multisample:
    input:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/genotypes/{sample}-multiallelic.vcf.gz', sample=samples)
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/genotypes/multisample-multiallelic.vcf')
    resources:
        mem_total_mb=180000,
        runtime_hrs=71,
    shell:
        '''
        set +u
        source ~/.bashrc
        conda activate base
        set -u
        bcftools merge --no-version -o {output} {input}
        '''


# get the biallelic multisample VCF
rule convert_multisample_vcf_biallelic:
    input:
        vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/genotypes/multisample-multiallelic.vcf',
        biallelic_vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/panel/giggles-ready_biallelic.vcf.gz'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/genotypes/multisample-biallelic.vcf')
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=3000
    shell:
        "cat {input.vcf} | python scripts/convert-to-biallelic.py {input.biallelic_vcf} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}"


### Quality Check ###

# calculate SV count
rule sv_count:
    input:
        metadata='/gpfs/project/projects/medbioinf/users/spani/files/other/1000GP/igsr_sample_data.tsv',
        vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/genotypes/multisample-biallelic.vcf.gz',
        script='scripts/plot-sv-count-boxplot.py'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/sv-counts-boxplot/het_count_population-wise.png',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/sv-counts-boxplot/het_count_sample-wise.png',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/sv-counts-boxplot/hom_count_population-wise.png',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/sv-counts-boxplot/hom_count_sample-wise.png',
        plot_data='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/sv-counts-boxplot/plot-sv-counts-boxplot.tsv'
    log:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/sv-counts-boxplot/plot-sv-counts-boxplot.log'
    params:
        out='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/sv-counts-boxplot/'
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=5000,
        runtime_hrs=24,
        runtime_min=1
    shell:
        'python {input.script} -meta {input.metadata} -vcf {input.vcf} -output {params.out} 2> {log}'

# calculate SV length and boxplot it
rule sv_length:
    input:
        metadata='/gpfs/project/projects/medbioinf/users/spani/files/other/1000GP/igsr_sample_data.tsv',
        vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/genotypes/multisample-biallelic.vcf.gz',
        script='scripts/plot-sv-length-boxplot.py'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/sv-length-boxplot/het_length_population-wise.png',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/sv-length-boxplot/het_length_sample-wise.png',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/sv-length-boxplot/hom_length_population-wise.png',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/sv-length-boxplot/hom_length_sample-wise.png',
        plot_data='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/sv-length-boxplot/plot-sv-length-boxplot.tsv'
    log:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/sv-length-boxplot/plot-sv-length-boxplot.log'
    params:
        out='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/sv-length-boxplot/'
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=5000,
        runtime_hrs=24,
        runtime_min=1
    shell:
        'python {input.script} -meta {input.metadata} -vcf {input.vcf} -output {params.out} 2> {log}'

# collect vcf stats
rule collect_vcf_stats:
    input:
        metadata='/gpfs/project/projects/medbioinf/users/spani/files/other/1000GP/igsr_sample_data.tsv',
        panel='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/panel/giggles-ready_biallelic.vcf.gz',
        callset='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/genotypes/multisample-biallelic.vcf.gz',
        script='scripts/collect-vcf-stats.py'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/hwe-af/callset-stats.tsv')
    log:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/hwe-af/callset-stats.log'
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=5000,
        runtime_hrs=71
    shell:
        'python {input.script} -meta {input.metadata} -panel {input.panel} -callset {input.callset} > {output} 2> {log}'

# adding bubble ids to the vcf stats table
rule add_bub_ids:
    input:
        table='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/hwe-af/callset-stats.tsv',
        multi_panel='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/panel/giggles-ready_multiallelic.vcf.gz',
        script='scripts/add-bub-info.py'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/hwe-af/variant-stats.tsv'
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=5000,
        runtime_hrs=3
    shell:
        'python {input.script} -table {input.table} -panel {input.multi_panel} -output {output}'

# plot af and hwe related statistics
rule calc_and_plot_statistics:
    input:
        table='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/hwe-af/variant-stats.tsv',
        script='scripts/plot-vcf-stats.py'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/hwe-af/hwe.png',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/hwe-af/af.png',
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/plots/hwe-af/hwe.{pop}.png', pop=['AFR', 'AMR', 'EAS', 'EUR', 'SAS']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/plots/hwe-af/af.{pop}.png', pop=['AFR', 'AMR', 'EAS', 'EUR', 'SAS']),
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/hwe-af/hwe.SVonly.png',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/hwe-af/af.SVonly.png',
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/plots/hwe-af/hwe.{pop}.SVonly.png', pop=['AFR', 'AMR', 'EAS', 'EUR', 'SAS']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/plots/hwe-af/af.{pop}.SVonly.png', pop=['AFR', 'AMR', 'EAS', 'EUR', 'SAS']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/plots/hwe-af/hwe.{vtype}.png', vtype=['INS', 'DEL', 'COMPLEX']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/plots/hwe-af/af.{vtype}.png', vtype=['INS', 'DEL', 'COMPLEX']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/plots/hwe-af/hwe.{vtype}.SVonly.png', vtype=['INS', 'DEL', 'COMPLEX']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/plots/hwe-af/af.{vtype}.SVonly.png', vtype=['INS', 'DEL', 'COMPLEX']),
        txt='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/hwe-af/output.stdout'
    params:
        outdir='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/hwe-af/'
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=5000,
        runtime_hrs=3
    shell:
        'python {input.script} -table {input.table} -output {params.outdir} > {output.txt}'

# plot sv discovery growth curve (audano curve), sv count distribution, and log(#svs) vs log(#ac) curve
rule plot_qc_curves:
    input:
        metadata='/gpfs/project/projects/medbioinf/users/spani/files/other/1000GP/igsr_sample_data.tsv',
        callset='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/genotypes/multisample-biallelic.vcf.gz',
        script='scripts/plot-qc-curves.py'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/qc/audano-curve.png',
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/plots/qc/average_count-size-distribution_{sv_type}_{ranges}.png', sv_type=['COMPLEX', 'DEL', 'INS'], ranges=['0-1kbp', '1kbp-10kbp', '10kbp-100kbp', '100kbp-1Mbp']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/plots/qc/total_count-size-distribution_{sv_type}_{ranges}.png', sv_type=['COMPLEX', 'DEL', 'INS'], ranges=['0-1kbp', '1kbp-10kbp', '10kbp-100kbp', '100kbp-1Mbp']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/plots/qc/stacked-average_count-size-distribution_{ranges}.png', ranges=['0-1kbp', '1kbp-10kbp', '10kbp-100kbp', '100kbp-1Mbp']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/plots/qc/stacked-total_count-size-distribution_{ranges}.png', ranges=['0-1kbp', '1kbp-10kbp', '10kbp-100kbp', '100kbp-1Mbp']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/plots/qc/rausch-curve_{v_set}-{v_type}.png', v_set = ['All', 'SV'], v_type = ['All', 'COMPLEX', 'DEL', 'INS'])
    params:
        out='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/plots/qc/'
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=5000
    shell:
        'python {input.script} -vcf {input.callset} -meta {input.metadata} -output {params.out}'

