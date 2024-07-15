# list of parameters with their values
par1 = []
par2 = []

# global wildcard constraints
wildcard_constraints:
    parameter1 = '|'.join(par1),
    parameter2 = '|'.join(par2)

# align chimp assemblies to rGFA using various parameters
rule align_assemblies_experiments:
    input:
        assembly=config['path_to_assembly'],
        graph=config['path_to_rgfa']
    output:
        'results/assembly-to-graph-alignment-experiments/alignment-{}-{}-{}.gaf'
    log:
        'results/assembly-to-graph-alignment-experiments/alignment.log'
    threads: 22
    resources:
        runtime_hrs=0,
        runtime_min=20,
        mem_total_mb=256
    conda:
        '../envs/alignment.yml'
    shell:
        'GraphAligner -g {input.graph} -f {input.assembly} -a {output} -t {threads} > {log}'

# get statistics using gaftools stats
rule alignment_statistics:
    input:
        'results/assembly-to-graph-alignment-experiments/alignment-{}-{}-{}.gaf'
    output:
        'results/assembly-to-graph-alignment-experiments/alignment-{}-{}-{}.stats'
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

# create fofn to collate statistics
rule fofn_statistics:
    input:
        expand('results/assembly-to-graph-alignment-experiments/alignment-{}-{}-{}.stats', ......parameters......)
    output:
        'results/assembly-to-graph-alignment-experiments/stats.fofn'
    run:
        from os.path import abspath
        f = open(output[0], 'w')
        for name in input:
            print(abspath(name), file=f)
        f.close()

# post-process all the stat files to get the relevant statistics together
rule collate_statistics:
    input:
        'results/assembly-to-graph-alignment-experiments/stats.fofn'
    output:
        'results/assembly-to-graph-alignment-experiments/collated-statistics.tsv'
    run:
        writer = open(output[0], 'w')
        reader = None
        with open(input[0], 'r') as fofn:
            for name in fofn:
                reader = open(name, 'r')
                for line in reader:
                    if line.startwith('Total aligned bases'):
                        print(name, line.split(' ')[-1], sep = '\t', file = writer)
                reader.close()
        writer.close()
