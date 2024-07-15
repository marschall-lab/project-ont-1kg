# align chimp assemblies to rGFA
rule align_assemblies:
    input:
        assembly=config['path_to_assembly'],
        graph=config['path_to_rgfa']
    output:
        'results/assembly-to-graph-alignment/alignment.gaf'
    log:
        'results/assembly-to-graph-alignment/alignment.log'
    threads: 22
    resources:
        runtime_hrs=,
        runtime_min=,
        mem_total_mb=
    conda:
        '../envs/alignment.yml'
    shell:
        'GraphAligner -g {input.graph} -f {input.assembly} -a {output} -t {threads} > {log}'

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
        graph=config['path_to_rgfa']
    output:
        'results/assembly-vcf/chimp.vcf'
    log:
        'results/assembly-vcf/chimp.log'
    resources:
        runtime_hrs=0,
        runtime_min=59,
        mem_total_mb=5000
    shell:
        '''
        set +u
        source ~/.bashrc
        conda activate giggles-env
        set -u
        giggles prepare-vcf --gfa {input.graph} --haploid {input.fofn} 1> {output} 2> {log}
        '''
    