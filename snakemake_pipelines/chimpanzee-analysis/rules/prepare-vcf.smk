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
        assembly=config['path_to_assembly'],
        graph=config['path_to_rgfa']
    output:
        'results/assembly-to-graph-alignment/alignment.gaf'
    log:
        'results/assembly-to-graph-alignment/alignment.log'
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
        'results/assembly-to-graph-alignment/alignment.gaf'
    output:
        'results/assembly-to-graph-alignment/alignment.stats'
    resources:
        runtime_hrs=0,
        runtime_min=20,
        mem_total_mb=256
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
        'results/assembly-to-graph-alignment/alignment.gaf'
    output:
        'results/assembly-to-graph-alignment/fofn.txt'
    run:
        from os.path import abspath
        f = open(output[0], 'w')
        for name in input:
            print(abspath(name), file=f)
        f.close()

# process the GAF using the rGFA to find the alleles
rule prepare_vcf:
    input:
        fofn='results/assembly-to-graph-alignment/fofn.txt',
        stats='results/assembly-to-graph-alignment/alignment.stats',
        graph=config['path_to_rgfa']
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
        giggles prepare_vcf --gfa {input.graph} --haploid {input.fofn} 1> {output} 2> {log}
        '''
    