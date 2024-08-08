# list of parameters with their values
seeds_mxm_widowsize = ['5000']
seeds_mxm_length = ['20']
seeds_mem_count = ['1000000']
bandwidth = ['15']
multimap_score_fraction = ['0.99']
precise_clipping = ['0.95']
min_alignment_score = ['5000']
clip_ambiguous_ends = ['100']
overlap_incompatible_cutoff = ['0.1']
max_trace_count = ['5']

# global wildcard constraints
wildcard_constraints:
    par1 = '|'.join(seeds_mxm_widowsize),
    par2 = '|'.join(seeds_mxm_length),
    par3 = '|'.join(seeds_mem_count),
    par4 = '|'.join(bandwidth),
    par5 = '|'.join(multimap_score_fraction),
    par6 = '|'.join(precise_clipping),
    par7 = '|'.join(min_alignment_score),
    par8 = '|'.join(clip_ambiguous_ends),
    par9 = '|'.join(overlap_incompatible_cutoff),
    par10 = '|'.join(max_trace_count)

# chop the references into synthetic long reads
rule extract_synthetic_reads:
    input:
        assembly=config['path_to_assembly']
    output:
        'results/chimp-long-reads.fa'
    log:
        'results/chimp-long-reads.log'
    resources:
        runtime_hrs=1,
        runtime_min=20,
        mem_total_mb=5*1024
    shell:
        'python scripts/chop-reference-to-reads.py -size 4000 -overlap 500 -ref {input} 1> {output} 2> {log}'
    

# align chimp assemblies to rGFA using various parameters
rule align_assemblies_experiments:
    input:
        reads='results/chimp-long-reads.fa',
        graph=config['path_to_rgfa']
    output:
        'results/assembly-to-graph-alignment-experiments/alignment-{par1}-{par2}-{par3}-{par4}-{par5}-{par6}-{par7}-{par8}-{par9}-{par10}.gaf'
    log:
        'results/assembly-to-graph-alignment-experiments/alignment-{par1}-{par2}-{par3}-{par4}-{par5}-{par6}-{par7}-{par8}-{par9}-{par10}.log'
    threads: 24
    resources:
        runtime_hrs=48,
        runtime_min=20,
        mem_total_mb=150*1024
    conda:
        '../envs/alignment.yml'
    shell:
        '''
        GraphAligner --seeds-mxm-windowsize {wildcards.par1} \
                     --seeds-mxm-length {wildcards.par2} \
                     --seeds-mem-count {wildcards.par3} \
                     --bandwidth {wildcards.par4} \
                     --multimap-score-fraction {wildcards.par5} \
                     --precise-clipping {wildcards.par6} \
                     --min-alignment-score {wildcards.par7} \
                     --clip-ambiguous-ends {wildcards.par8} \
                     --overlap-incompatible-cutoff {wildcards.par9} \
                     --max-trace-count {wildcards.par10} \
                     -g {input.graph} \
                     -f {input.reads} \
                     -a {output} \
                     -t {threads} > {log}
        '''

# get statistics using gaftools stats
rule alignment_statistics:
    input:
        'results/assembly-to-graph-alignment-experiments/alignment-{par1}-{par2}-{par3}-{par4}-{par5}-{par6}-{par7}-{par8}-{par9}-{par10}.gaf'
    output:
        'results/assembly-to-graph-alignment-experiments/alignment-{par1}-{par2}-{par3}-{par4}-{par5}-{par6}-{par7}-{par8}-{par9}-{par10}.stats'
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
        expand('results/assembly-to-graph-alignment-experiments/alignment-{par1}-{par2}-{par3}-{par4}-{par5}-{par6}-{par7}-{par8}-{par9}-{par10}.stats',
                    par1 = seeds_mxm_widowsize,
                    par2 = seeds_mxm_length,
                    par3 = seeds_mem_count,
                    par4 = bandwidth,
                    par5 = multimap_score_fraction,
                    par6 = precise_clipping,
                    par7 = min_alignment_score,
                    par8 = clip_ambiguous_ends,
                    par9 = overlap_incompatible_cutoff,
                    par10 = max_trace_count)
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
        header = ['seeds_mxm_widowsize',
                    'seeds_mxm_length',
                    'seeds_mem_count',
                    'bandwidth',
                    'multimap_score_fraction',
                    'precise_clipping',
                    'min_alignment_score',
                    'clip_ambiguous_ends',
                    'overlap_incompatible_cutoff',
                    'max_trace_count']
        print('\t'.join(header)+'\tstat', file=writer)
        with open(input[0], 'r') as fofn:
            for name in fofn:
                reader = open(name, 'r')
                for line in reader:
                    if line.startwith('Total aligned bases'):
                        print('\t'.join(name.split('/')[-1].split('.')[0].split('-')[1:]), line.split(' ')[-1], sep = '\t', file = writer)
                reader.close()
        writer.close()
