# extract fasta from cram
rule cram_to_fasta:
    input:
        cram='/gpfs/project/projects/medbioinf/data/share/globus/1000g-ont/hg38/{sample}.hg38.cram',
        ref='/gpfs/project/projects/medbioinf/users/spani/files/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta'
    output:
        fasta='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/data/fasta/{sample}.fasta.gz',
        fai='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/data/fasta/{sample}.fasta.gz.fai',
        gzi='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/data/fasta/{sample}.fasta.gz.gzi'
    conda:
        "../envs/basic.yml"
    resources:
        runtime_hrs=20,
        mem_total_mb=20000
    priority: 1
    shell:
        '''
        seq_cache_populate.pl -root /gpfs/project/projects/medbioinf/users/spani/files/ref {input.ref}
        export REF_PATH=/gpfs/project/projects/medbioinf/users/spani/files/ref/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
        export REF_CACHE=/gpfs/project/projects/medbioinf/users/spani/files/ref/%2s/%2s/%s
        samtools fasta {input.cram} | bgzip -c > {output.fasta}
        samtools faidx {output.fasta}
        '''

# run minigraph alignment
rule minigraph_alignment:
    input:
        fasta='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/data/fasta/{sample}.fasta.gz',
        fai='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/data/fasta/{sample}.fasta.gz.fai',
        gzi='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/data/fasta/{sample}.fasta.gz.gzi',
        ref='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/ont-alignments/{sample}.gaf')
    log:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/ont-alignments/{sample}.alignment.log'
    resources:
        runtime_hrs=36,
        runtime_min=0,
        mem_total_mb=96000
    threads: 8
    shell:
        '/gpfs/project/projects/medbioinf/users/spani/packages/minigraph/minigraph --vc -cx lr {input.ref} {input.fasta} -t {threads} > {output}'


# gaf sorting
rule gaftools_sort:
    input:
        gaf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/ont-alignments/{sample}.gaf',
        gfa='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa'
    output:
        sorted_gaf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/ont-alignments/{sample}.sorted.gaf.gz',
        index='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/ont-alignments/{sample}.sorted.gaf.gz.gai'
    log:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{callset}/ont-alignments/{sample}.sorting.log'
    resources:
        runtime_hrs=54,
        runtime_min=0,
        mem_total_mb=20000
    shell:
        '''
        set +u
        source ~/.bashrc
        conda activate gaftools-dev
        set -u
        gaftools sort --bgzip --outgaf {output.sorted_gaf} {input.gaf} {input.gfa} 2> {log}
        '''