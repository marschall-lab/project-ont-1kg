configfile: 'config/config.yaml' 

# tag rGFA chromosome-wise
rule call_rGFA_bubbles:
    input:
        ref=lambda wildcards: config['callsets'][wildcards.callset]['gfa']
    output:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/bubble_calling_and_tagging/{{callset}}-{chr}.gfa', chr=chromosomes),
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/bubble_calling_and_tagging/{{callset}}-{chr}.csv', chr=chromosomes)
    params:
        out_dir='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/bubble_calling_and_tagging'
    resources:
        runtime_hrs=0,
        runtime_min=30,
        mem_total_mb=lambda wildcards, attempt: 10*1024 * attempt
    shell:
        '''
        set +u
        source ~/.bashrc
        conda activate gaftools-dev
        set -u
        gaftools order_gfa --with-sequence --outdir {params.out_dir} {input.ref}
        '''

# concat chromsome-wise tagged GFA
rule concat_tagged_GFA:
    input:
       lambda wildcards: expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{{callset}}/bubble_calling_and_tagging/%s-{chr}.gfa'%(config['callsets'][wildcards.callset]['gfa'].split("/")[-1][:-4]), chr=chromosomes)
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa'
    shell:
        'cat {input} > {output}'