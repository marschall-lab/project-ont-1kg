#extract gaf (to make gaftools sort faster(?))
rule decompress_gaf:
    input:
        '/gpfs/project/projects/medbioinf/data/share/globus/1000g-ont/gaf/{sample}.gaf.gz'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/gaf/{sample}.gaf')
    resources:
        runtime_hrs=2,
        mem_total_mb=20000
    shell:
        'gzip -d -c {input} > {output}'

#compress vcf
rule compress_vcf:
    input:
        "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/{filename}.vcf"
    output:
        vcf="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/{filename}.vcf.gz",
        tbi="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/{filename}.vcf.gz.tbi"
    shell:
        """
        bgzip -c {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# extract fasta from cram
rule cram_to_fasta:
    input:
        cram='/gpfs/project/projects/medbioinf/data/share/globus/1000g-ont/hg38/{sample}.hg38.cram',
        ref='/gpfs/project/projects/medbioinf/users/spani/files/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta'
    output:
        fasta="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz",
        fai="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz.fai",
        gzi="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz.gzi"
    conda:
        "../envs/basic.yml"
    resources:
        runtime_hrs=10,
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

# gaf sorting
rule gaftools_sort:
    input:
        gaf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/gaf/{sample}.gaf',
        gfa='/gpfs/project/projects/medbioinf/users/spani/files/gfa/HengLi/chm13-90c.r518_tagged.gfa'
    output:
        sorted_gaf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/gaf/{sample}.sorted.gaf.gz',
        index='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/gaf/{sample}.sorted.gaf.gz.gai'
    log:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/1000GP/svarp-giggles/chm13-90c.r518/data/gaf/{sample}.log'
    resources:
        runtime_hrs=24,
        runtime_min=0,
        mem_total_mb=10000
    shell:
        '''
        set +u
        source ~/.bashrc
        conda activate gaftools-dev
        set -u
        gaftools sort --bgzip --outgaf {output.sorted_gaf} {input.gaf} {input.gfa}
        '''

# prepare vcf panel
rule prepare_vcf:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/annotating-paths/vcf/chm13-90c.r518_filtered.vcf.gz'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel.vcf'
    conda:
        "../envs/basic.yml"
    resources:
        mem_total_mb=20000,
        runtime_hrs=0,
        runtime_min=30
    shell:
        "bcftools view --trim-alt-alleles -c1 {input} > {output}"