#compress vcf
rule compress_vcf:
    input:
        "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/{filename}.vcf"
    output:
        vcf = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/{filename}.vcf.gz",
        tbi = "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/{filename}.vcf.gz.tbi"
    priority: 1
    shell:
        """
        bgzip -c {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# extract fasta from cram
rule cram_to_fasta:
    input:
        cram = '/gpfs/project/projects/medbioinf/data/share/globus/1000g-ont/GRCh38/{sample}/alignments/{sample}.cram',
        ref= '/gpfs/project/projects/medbioinf/users/spani/files/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta'
    output:
        fasta = temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz"),
        fai = temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz.fai"),
        gzi = temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz.gzi")
    conda:
        "../envs/basic.yml"
    resources:
        runtime_hrs=10,
        mem_total_mb=20000
    shell:
        '''
        seq_cache_populate.pl -root /gpfs/project/projects/medbioinf/users/spani/files/ref {input.ref}
        export REF_PATH=/gpfs/project/projects/medbioinf/users/spani/files/ref/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
        export REF_CACHE=/gpfs/project/projects/medbioinf/users/spani/files/ref/%2s/%2s/%s
        samtools fasta {input.cram} | bgzip -c > {output.fasta}
        samtools faidx {output.fasta}
        '''

# gaf sorting
rule gaftools_sort:
    input:
        gaf = '/gpfs/project/projects/medbioinf/data/share/globus/1000g-ont/gaf/{sample}.gaf.gz',
        gfa = '/gpfs/project/projects/medbioinf/users/spani/files/gfa/HengLi/chm13-90c.r518_tagged.gfa'
    output:
        sorted_gaf = '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/gaf/{sample}.sorted.gaf.gz',
        index = '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/gaf/{sample}.sorted.gaf.gz.gai'
    log:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/1000GP/svarp-giggles/chm13-90c.r518/data/gaf/{sample}.log'
    resources:
        runtime_hrs=6,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 20000 * attempt
    shell:
        '''
        set +u
        source ~/.bashrc
        conda activate gaftools-dev
        set -u
        gaftools sort --bgzip --outgaf {output.sorted_gaf} {input.gaf} {input.gfa}
        '''