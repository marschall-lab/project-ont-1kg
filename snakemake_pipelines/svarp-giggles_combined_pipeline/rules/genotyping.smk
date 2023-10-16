include: './get-sample-list.smk'

# run genotyping
rule giggles:
    input:
        reads='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz',
        reads_gzi='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz.gzi',
        reads_fai='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz.fai',
        haplotag='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/haplotagging-results/GRCh38/{sample}/{sample}.tsv',
        alignment='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/gaf/{sample}.sorted.gaf.gz',
        alignment_index='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/gaf/{sample}.sorted.gaf.gz.gai',
        gfa='/gpfs/project/projects/medbioinf/users/spani/files/gfa/HengLi/chm13-90c.r518_tagged.gfa',
        vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel.vcf.gz'
    output:
        vcf=temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/{sample}.vcf")
    log:
        stderr="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/{sample}.stderr",
        stdout="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/{sample}.stdout"
    resources:
        mem_total_mb=50000,
        runtime_hrs=48,
        runtime_min=1
    priority: 3
    shell:
        """
        set +u
        source ~/.bashrc
        conda activate giggles-dev
        set -u
        giggles genotype --read-fasta {input.reads} --sample {wildcards.sample} --realign-mode edit -o {output.vcf} --rgfa {input.gfa} --haplotag-tsv {input.haplotag} {input.vcf} {input.alignment} > {log.stdout} 2> {log.stderr}
        """

# add allele decomposition info to the panel VCF and get the biallelic and multiallelic panel VCF.
rule annotate_panel:
    input:
        vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel.vcf.gz',
        gfa='/gpfs/project/projects/medbioinf/users/spani/files/gfa/HengLi/chm13-90c.r518_tagged.gfa'
    output:
        vcf_temp=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel-unzipped-tmp.vcf'),
        multi='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel-multiallelic.vcf',
        multi_tmp=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel-tmp.vcf'),
        biallelic='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel-biallelic.vcf',
        bi_tmp=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel-tmp_biallelic.vcf')
    log:
        "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel-annotate.log"
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=20000,
        runtime_hrs=5,
        runtime_min=59
    params:
        outname='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel-tmp'
    shell:
        """
        gzip -d -c {input.vcf} > {output.vcf_temp}
        python3 scripts/annotate_vcf.py -vcf {output.vcf_temp} -gfa {input.gfa} -o {params.outname} &> {log}
        cat {output.multi_tmp} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' > {output.multi}
        cat {output.bi_tmp} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' > {output.biallelic}
        """

# add allele decomposition info to the genotyped VCFs
rule add_allele_decomposition_info:
    input:
        sample_vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/{sample}.vcf.gz',
        panel_vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel-multiallelic.vcf'
    output:
        sample_vcf_temp=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/{sample}-tmp.vcf'),
        output_vcf=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/{sample}-multiallelic.vcf')
    resources:
        mem_total_mb=3000
    shell:
        '''
        gzip -d -c {input.sample_vcf} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' > {output.sample_vcf_temp}
        awk -F"\t" -f scripts/copy_columns.awk {input.panel_vcf} {output.sample_vcf_temp} > {output.output_vcf}
        '''
    

# get the biallelic and multiallelic VCFs for the genotypes
rule convert_sample_vcf_biallelic:
    input:
        sample_vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/{sample}-multiallelic.vcf',
        biallelic_vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel-biallelic.vcf.gz'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/{sample}-biallelic.vcf')
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=3000
    shell:
        "cat {input.sample_vcf} | python scripts/convert-to-biallelic.py {input.biallelic_vcf} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}"


# merge multiallelic records to give multisample vcf
rule merge_vcf_to_multisample:
    input:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/{sample}-multiallelic.vcf.gz', sample=samples)
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/multisample-multiallelic.vcf')
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=30000,
        runtime_hrs=24,
    shell:
        'bcftools merge --no-version -o {output} {input}'


# get the biallelic multisample VCF
rule convert_multisample_vcf_biallelic:
    input:
        sample_vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/multisample-multiallelic.vcf',
        biallelic_vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel-biallelic.vcf.gz'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/multisample-biallelic.vcf')
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=3000
    shell:
        "cat {input.sample_vcf} | python scripts/convert-to-biallelic.py {input.biallelic_vcf} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}"


# calculate SV count
rule sv_count:
    input:
        metadata='/gpfs/project/projects/medbioinf/users/spani/files/other/1000GP/igsr_sample_data.tsv',
        vcf=expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/{sample}-biallelic.vcf',sample=samples)
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/het_count_population-wise.png',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/het_count_sample-wise.png',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/hom_count_population-wise.png',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/hom_count_sample-wise.png',
        plot_data='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/sv_count.tsv'
    log:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/sv_count.log'
    params:
        inp_vcfs=lambda wildcards, input: ','.join(input.vcf),
        out='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/'
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=5000,
        runtime_hrs=24,
        runtime_min=1
    shell:
        'python scripts/sv_count.py -meta {input.metadata} -vcf {params.inp_vcfs} -output {params.out} 2> {log}'

# collect vcf stats
rule collect_vcf_stats:
    input:
        metadata='/gpfs/project/projects/medbioinf/users/spani/files/other/1000GP/igsr_sample_data.tsv',
        panel='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel-biallelic.vcf.gz',
        callset='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/multisample-biallelic.vcf.gz'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/callset-stats.tsv')
    log:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/callset-stats.log'
    conda:
        '../envs/cyvcf2.yml'
    resources:
        mem_total_mb=5000,
        runtime_hrs=18
    shell:
        'python scripts/collect-vcf-stats.py -meta {input.metadata} -panel {input.panel} -callset {input.callset} > {output} 2> {log}'

# adding bubble ids to the vcf stats table
rule add_bub_ids:
    input:
        table='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/callset-stats.tsv',
        multi_panel='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel-multiallelic.vcf'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/variant-stats.tsv'
    conda:
        '../envs/cyvcf2.yml'
    resources:
        mem_total_mb=5000,
        runtime_hrs=3
    shell:
        'python scripts/add-bub-info.py -table {input.table} -panel {input.multi_panel} -output {output}'

# plot statistics
rule plot_statistics:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/variant-stats.tsv'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/hwe.png',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/af.png',
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/hwe.{pop}.png', pop=['AFR', 'AMR', 'EAS', 'EUR', 'SAS']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/af.{pop}.png', pop=['AFR', 'AMR', 'EAS', 'EUR', 'SAS']),
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/hwe.SVonly.png',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/af.SVonly.png',
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/hwe.{pop}.SVonly.png', pop=['AFR', 'AMR', 'EAS', 'EUR', 'SAS']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/af.{pop}.SVonly.png', pop=['AFR', 'AMR', 'EAS', 'EUR', 'SAS']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/hwe.{vtype}.png', vtype=['INS', 'DEL', 'COMPLEX']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/af.{vtype}.png', vtype=['INS', 'DEL', 'COMPLEX']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/hwe.{vtype}.SVonly.png', vtype=['INS', 'DEL', 'COMPLEX']),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/af.{vtype}.SVonly.png', vtype=['INS', 'DEL', 'COMPLEX']),
    params:
        outdir='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/plots/'
    conda:
        '../envs/cyvcf2.yml'
    shell:
        'python scripts/plot-vcf-stats.py -table {input} -output {params.outdir}'