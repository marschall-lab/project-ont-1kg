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

# extract bubbles from gfa
rule get_bubbles_from_gfa:
    input:
        'results/graph-tagged.gfa'
    output:
        'results/gfa-bubbles.bed'
    shell:
        'python scripts/extract-bubbles.py -gfa {input} > {output}'

# evaluate recovery of SVs
rule evaluate_vcf_recovery:
    input:
        bubbles_bed='results/gfa-bubbles.bed',
        panel_vcf=config['path_to_panel_vcf']
        callset_vcf=config['path_to_callset_vcf']
    output:
        'results/recovery-stats.txt'
    params:
        'results/'
    shell:
        'python scripts/evaluate-sv-recovery.py -bed {input.bubbles_bed} -panel {input.panel_vcf} -callset {input.callset_vcf} -outdir {params}'