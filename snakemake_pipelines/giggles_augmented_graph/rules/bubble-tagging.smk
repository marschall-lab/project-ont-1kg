# tag rGFA chromosome-wise
rule call_rGFA_bubbles:
    input:
        ref=lambda wildcards: config['callsets'][wildcards.callset]['gfa']
    output:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{{callset}}/bubble_calling_and_tagging/{chr}.gfa', chr=chromosomes),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{{callset}}/bubble_calling_and_tagging/{chr}.csv', chr=chromosomes)
    params:
        out_dir='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{{callset}}/bubble_calling_and_tagging'
    resources:
        runtime_hrs=0,
        runtime_min=30,
        mem_total_mb=lambda wildcards, attempt: 10*1024 * attempt
    shell:
        '''
        module load Python/3.11.4
        gaftools order_gfa --with-sequence --outdir {params.out_dir} {input.ref}
        module unload Python/3.11.4
        '''

# concat chromsome-wise tagged GFA
rule concat_tagged_GFA:
    input:
       expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{{callset}}/bubble_calling_and_tagging/{chr}.gfa', chr=chromosomes)
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa'
    shell:
        'cat {input} > {output}'