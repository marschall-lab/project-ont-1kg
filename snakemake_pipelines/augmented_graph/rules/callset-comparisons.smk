
include: './get-sample-list.smk'

if config['pilot']:
    samples.sort()
    samples=samples[0:2]

max_af = 1
min_af = 0
chromosomes='chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX'.split(',')
# add nygc here when the liftover vcf is ready
sources='pangenie,giggles,pangenie_panel,giggles_panel'.split(',')

wildcard_constraints:
    chr='|'.join(chromosomes),
    sample='|'.join(samples),
    source='|'.join(sources)

### extract sample from nygc phased panel ###
# extract sample from chromosome-wise vcf
rule nygc_extract_sample_chr_vcf:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/files/vcf/1000GP/NYGC/phased/1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/nygc-{sample}-{chr}.vcf')
    conda:
        '../envs/basic.yml'
    shell:
        'bcftools view --samples {wildcards.sample} {input} > {output}'

# join chromosome-wise vcfs for a sample
rule nygc_merge_chr_vcfs:
    input:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/callset-comparison/vcfs/nygc-{{sample}}-{chr}.vcf', chr=chromosomes)
    output:
        vcf=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/nygc-{sample}.vcf.gz'),
        index=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/nygc-{sample}.vcf.gz.tbi')
    conda:
        '../envs/basic.yml'
    shell:
        '''
        bcftools concat -o {output.vcf} -Oz {input}
        tabix -p vcf {output.vcf}
        '''

# extract sample from pangenie genotyped vcf
rule pangenie_create_sample_vcf:
    input:
        '/gpfs/project/projects/medbioinf/users/ebler/long-read-1kg/pangenie-genotypes/results/pangenie-{sample}_genotyping_biallelic.vcf.gz'
    output:
        vcf=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/pangenie-{sample}.vcf.gz'),
        index=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/pangenie-{sample}.vcf.gz.tbi')
    conda:
        '../envs/basic.yml'
    shell:
        '''
        bcftools view {input} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        '''

# extract sample from giggles genotyped vcf
rule giggles_create_sample_vcf:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/genotypes/{sample}-biallelic.vcf.gz'
    output:
        vcf=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/giggles-{sample}.vcf.gz'),
        index=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/giggles-{sample}.vcf.gz.tbi')
    conda:
        '../envs/basic.yml'
    shell:
        '''
        bcftools view {input} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        '''

# extract sample vcf from pangenie panel (only for HG01258)
rule pangenie_panel_extract_sample:
    input:
        '/gpfs/project/projects/medbioinf/users/ebler/long-read-1kg/pangenie-genotypes/input/chm13_cactus_filtered_ids_biallelic.vcf'
    output:
        vcf=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/pangenie_panel-HG01258.vcf.gz'),
        index=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/pangenie_panel-HG01258.vcf.gz.tbi')
    conda:
        '../envs/basic.yml'
    shell:
        '''
        bcftools view --samples HG01258 {input} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        '''

# extract vcf from giggles panel (only for HG01258)
rule giggles_panel_extract_sample:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/panel/giggles-ready_biallelic.vcf.gz'
    output:
        vcf=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/giggles_panel-HG01258.vcf.gz'),
        index=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/giggles_panel-HG01258.vcf.gz.tbi')
    conda:
        '../envs/basic.yml'
    shell:
        '''
        bcftools view --samples HG01258 {input} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        '''

# basic counting stats of the vcfs
rule vcf_stats:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/{source}-{sample}.vcf.gz'
    output:
        "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/{source}-{sample}.stats"
    conda:
        "../envs/basic.yml"
    shell:
        "zcat {input} | python scripts/set-pass.py | bcftools view -f PASS --min-af {min_af} --max-af {max_af} | python scripts/vcf_stats.py > {output}"

# add tag of variant types
rule add_tags:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/{source}-{sample}.vcf.gz'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/{source}-{sample}-tagged.vcf')
    conda:
        "../envs/basic.yml"
    shell:
        "zcat {input} | python scripts/set-pass.py | bcftools view -f PASS --min-af {min_af} --max-af {max_af} | python3 scripts/add-svtags.py > {output}"


# intersecting the vcf files and creating upset plots for HG01258 (sample in graph)
rule intersect_HG01258_vcfs:
    input:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/callset-comparison/vcfs/{source}-HG01258-tagged.vcf', source=sources)
    output:
        tsv="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/in-graph/HG01258-intersection.tsv",
        vcf="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/in-graph/HG01258-intersection.vcf",
        pdf="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/in-graph/HG01258-intersection.pdf",
        plot="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/in-graph/HG01258-intersection-upsetplot.pdf"
    conda:
        "../envs/basic.yml"
    log:
        intersect="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/in-graph/HG01258-intersection.log",
        plot="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/in-graph/HG01258-plotting.log"
    params:
        names=sources,
        columns=["in_" + s for s in sources]
    shell:
        '''
        python scripts/intersect_callsets.py intersect -c {input} -n {params.names} -t {output.tsv} -v {output.vcf} -p {output.pdf} &> {log.intersect}
        python scripts/plot-comparison-upset.py -t {output.tsv} -o {output.plot} -n {params.columns} &> {log.plot}
        '''

# intersecting the vcf files and creating upset plots for samples not in the graph
rule intersect_outsample_vcfs:
    input:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/callset-comparison/vcfs/{source}-{{sample}}-tagged.vcf', source=['pangenie', 'giggles'])
    output:
        tsv="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/out-graph/{sample}-intersection.tsv",
        vcf="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/out-graph/{sample}-intersection.vcf",
        pdf="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/out-graph/{sample}-intersection.pdf",
        plot="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/out-graph/{sample}-intersection-upsetplot.pdf"
    conda:
        "../envs/basic.yml"
    log:
        intersect="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/out-graph/{sample}-intersection.log",
        plot="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/out-graph/{sample}-plotting.log"
    params:
        names=['nygc', 'pangenie', 'giggles'],
        columns=["in_" + s for s in ['nygc', 'pangenie', 'giggles']]
    shell:
        '''
        python scripts/intersect_callsets.py intersect -c {input} -n {params.names} -t {output.tsv} -v {output.vcf} -p {output.pdf} &> {log.intersect}
        python scripts/plot-comparison-upset.py -t {output.tsv} -o {output.plot} -n {params.columns} &> {log.plot}
        '''
'''
### truvari comparisons ###

# compare giggles genotypes to pangenie genotypes
rule truvari_pangenie_genotypes_compare:
    input:
        giggles='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/giggles-{sample}.vcf.gz',
        pangenie='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/callset-comparison/vcfs/pangenie-{sample}.vcf.gz'
    output:




# compare giggles genotypes to pangenie panel (only for HG01258)



# compare giggles genotypes to nygc panel 
'''