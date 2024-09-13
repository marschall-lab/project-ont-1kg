# align chimp assemblies to rGFA using GraphAligner
# rule graphaligner_align_assemblies:
#     input:
#         assembly=config['path_to_assembly'],
#         graph=config['path_to_rgfa']
#     output:
#         'results/assembly-to-graph-alignment/alignment.gaf'
#     log:
#         'results/assembly-to-graph-alignment/alignment.log'
#     threads: 24
#     resources:
#         runtime_hrs=48,
#         runtime_min=59,
#         mem_total_mb=150*1024
#     conda:
#         '../envs/alignment.yml'
#     shell:
#         'GraphAligner -g {input.graph} -f {input.assembly} -a {output} -t {threads} > {log}'

# align chimp reference to rGFA using minigraph
rule minigraph_align_assemblies:
    input:
        assembly='results/chimp-long-reads.fa',
        graph=config['path_to_rgfa']
    output:
        'results/assembly-to-graph-alignment/chimp.gaf'
    log:
        'results/assembly-to-graph-alignment/chimp.log'
    threads: 24
    resources:
        runtime_hrs=24,
        runtime_min=59,
        mem_total_mb=80*1024
    conda:
        '../envs/alignment.yml'
    shell:
        'minigraph --vc -cx lr {input.graph} {input.assembly} -t {threads} > {output} 2> {log}'

# get statistics using gaftools stats
rule alignment_statistics:
    input:
        'results/assembly-to-graph-alignment/chimp.gaf'
    output:
        'results/assembly-to-graph-alignment/chimp.stats'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=2048
    shell:
        '''
        set +u
        source ~/.bashrc
        conda activate gaftools-env
        set -u
        gaftools stat --cigar {input} > {output}
        '''

# create FOFN for the prepare vcf script
rule make_list_of_file_names:
    input:
        'results/assembly-to-graph-alignment/chimp.gaf'
    output:
        'results/assembly-to-graph-alignment/fofn.txt'
    run:
        from os.path import abspath
        f = open(output[0], 'w')
        for name in input:
            print(abspath(name), file=f)
        f.close()

chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']

# tag rGFA chromosome-wise
rule call_rGFA_bubbles:
    input:
        ref=config['path_to_rgfa']
    output:
        expand('results/rgfa-tagging/%s-{chr}.gfa'%(config['path_to_rgfa'].split("/")[-1][:-4]), chr=chromosomes),
        'results/rgfa-tagging/%s-complete.gfa'%(config['path_to_rgfa'].split("/")[-1][:-4]),
        expand('results/rgfa-tagging/%s-{chr}.csv'%(config['path_to_rgfa'].split("/")[-1][:-4]), chr=chromosomes)
    params:
        out_dir='results/rgfa-tagging'
    resources:
        runtime_hrs=0,
        runtime_min=30,
        mem_total_mb=lambda wildcards, attempt: 10*1024 * attempt
    shell:
        '''
        set +u
        source ~/.bashrc
        conda activate gaftools-env
        set -u
        gaftools order_gfa --with-sequence --outdir {params.out_dir} {input.ref}
        '''

# concat chromsome-wise tagged GFA
rule concat_tagged_GFA:
    input:
        expand('results/rgfa-tagging/%s-{chr}.gfa'%(config['path_to_rgfa'].split("/")[-1][:-4]), chr=chromosomes)
    output:
        'results/graph-tagged.gfa'
    shell:
        'cat {input} > {output}'    

# process the GAF using the rGFA to find the alleles
rule prepare_vcf:
    input:
        fofn='results/assembly-to-graph-alignment/fofn.txt',
        stats='results/assembly-to-graph-alignment/chimp.stats',
        graph='results/rgfa-tagging/%s-complete.gfa'%(config['path_to_rgfa'].split("/")[-1][:-4])
    output:
        'results/assembly-vcf/chimp.vcf'
    log:
        'results/assembly-vcf/chimp.log'
    resources:
        runtime_hrs=0,
        runtime_min=59,
        mem_total_mb=5*1024
    shell:
        '''
        set +u
        source ~/.bashrc
        conda activate giggles-env
        set -u
        giggles prepare_vcf --gfa {input.graph} --haploid {input.fofn} --keep-all-records 1> {output} 2> {log}
        '''

# create a BED file with the ancestral allele annotations
rule annotate_ancestral_allele:
    input:
        'results/assembly-vcf/chimp.vcf'
    output:
        'results/ancestral-allele-annotations.bed'
    log:
        'results/ancestral-allele-annotations.log'
    shell:
        'python scripts/annotate-ancestral-allele.py -vcf {input} 1> {output} 2> {log}'

# create map between sv alleles and the bubble ids
rule bub_sv_map:
    input:
        config['path_to_multiallelic_callset']
    output:
        'results/bub-sv-map.tsv'
    conda:
        '../envs/analysis.yml'
    shell:
        'python scripts/map-bub-ids-to-allele-ids.py -vcf {input} > {output}'

# unzip phased callset vcf
rule unzip_phased_vcf:
    input:
        config['path_to_callset_vcf']
    output:
        temp('results/tmp/phased-callset.vcf')
    shell:
        'gzip -d -c {input} > {output}'

# match the ancestral alleles to the sv alleles and modify the vcf
rule match_ancestral_allele:
    input:
        bed='results/ancestral-allele-annotations.bed',
        bmap='results/bub-sv-map.tsv',
        vcf='results/tmp/phased-callset.vcf'
    output:
        'results/annotated-callset.vcf'
    log:
        'results/annotated-callset.log'
    shell:
        'python scripts/match-ancestral-alleles.py -bed {input.bed} -map {input.bmap} -vcf {input.vcf} 1> {output} 2> {log}'