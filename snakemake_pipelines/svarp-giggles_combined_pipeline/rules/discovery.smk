# run svarp
rule svarp:
    input:
        reads='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz',
        reads_gzi='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz.gzi',
        reads_fai='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz.fai',
        haplotag='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/haplotagging-results/GRCh38/{sample}/{sample}.tsv',
        alignment='/gpfs/project/projects/medbioinf/data/share/globus/1000g-ont/gaf/{sample}.gaf.gz',
        gfa='/gpfs/project/projects/medbioinf/users/spani/files/gfa/1000GP/chm13-90c.r518_tagged.gfa'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_H1.fa',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_H2.fa',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_untagged.fa'
    params:
        outdir='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}'
    log:
        stdout='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/log/{sample}.stdout',
        stderr='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/log/{sample}.stderr'
    conda:
        "../envs/svarp.yml"
    resources:
        mem_total_mb=1024*400,
        runtime_hrs=24*5,
        runtime_min=1
    priority: 2
    shell:
        '''
        module load Python/3.11.3
        module load Minimap2/2.17
        module load SamTools/1.6
        module load gcc/10.2.0
        export PATH=$PATH:/gpfs/project/projects/medbioinf/users/spani/packages/wtdbg2
        /gpfs/project/projects/medbioinf/projects/1000g-ont/svarp/build/svarp -a {input.alignment} -g {input.gfa} --fasta {input.reads} -i {wildcards.sample} --phase {input.haplotag} -o {params.outdir} 2> {log.stderr} 1> {log.stdout}
        module unload gcc/10.2.0
        module unload Minimap2/2.17
        module unload Python/3.11.3
        module unload SamTools/1.6
        '''

# map svtigs with minimap2
rule svtigs_minimap2:
    input:
        ref='/gpfs/project/projects/medbioinf/users/ebler/long-read-1kg/data/1KG_ONT_VIENNA_hg38.fa',
        fasta='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_{haplotype}.fa'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_{haplotype}.sam')
    wildcard_constraints:
        haplotype='H1|H2'
    conda:
        '../envs/svarp_processing.yml'
    resources:
        mem_total_mb=30000
    threads: 2
    priority:3
    shell:
        'minimap2 -a -x asm5 --cs -r2k -t {threads} {input.ref} {input.fasta}  > {output}'

# sort aligned SAM file
rule sort_svtig_alignment:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_{haplotype}.sam'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_{haplotype}.sorted.bam')
    wildcard_constraints:
        haplotype='H1|H2'
    conda:
        '../envs/svarp_processing.yml'
    threads: 2
    priority:3
    shell:
        'samtools sort -m4G -@{threads} -o {output} {input}'

# index sorted BAM
rule index_sorted_svtig_alignment:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_{haplotype}.sorted.bam'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_{haplotype}.sorted.bam.bai')
    wildcard_constraints:
        haplotype='H1|H2'
    conda:
        '../envs/svarp_processing.yml'
    shell:
        'samtools index {input}'

# run svim-asm to get the vcf
rule run_svim_asm:
    input:
        ref='/gpfs/project/projects/medbioinf/users/ebler/long-read-1kg/data/1KG_ONT_VIENNA_hg38.fa',
        bam1='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_H1.sorted.bam',
        bam1_index='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_H1.sorted.bam.bai',
        bam2='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_H2.sorted.bam',
        bam2_index='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_H2.sorted.bam.bai'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/variants.vcf',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/sv-lengths.png'
    params:
        outdir='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/'
    conda:
        '../envs/svarp_processing.yml'
    priority: 3
    shell:
        'svim-asm diploid --sample {wildcards.sample} {params.outdir} {input.bam1} {input.bam2} {input.ref}'

rule pav_pipeline:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_H1.fa',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_H2.fa'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav_hg38/run.complete',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav_hg38/pav_svtigs.vcf.gz',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav_t2t/run.complete',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav_t2t/pav_svtigs.vcf.gz'
    params:
        dir='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/'
    threads: 8
    resources:
        runtime_hrs=12,
        mem_total_mb=70000
    shell:
        '''
        module load Snakemake/7.8.5
        module load Singularity
        cd {params.dir}
        mkdir -p pav_hg38
        mkdir -p pav_t2t
        snakemake -s /gpfs/project/projects/medbioinf/users/spani/scripts/1000g-ont/ont-1kg/snakemake_pipelines/svarp-giggles_combined_pipeline/rules/pav.smk -c {threads} --use-singularity --config sample={wildcards.sample} -F
        rm -r pav_hg38/data
        rm -r pav_hg38/temp
        rm -r pav_hg38/results
        rm -r pav_t2t/data
        rm -r pav_t2t/temp
        rm -r pav_t2t/results
        '''

rule process_pav_output:
    input:
        ref='/gpfs/project/projects/medbioinf/users/spani/files/ref/1KG_ONT_VIENNA_{ref}.fa',
        vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav_{ref}/pav_svtigs.vcf.gz'
    output:
        final='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav_{ref}/pav_svtigs.vcf.gz_merged.vcf'
    conda:
        '../envs/svarp_processing.yml'
    wildcard_constraints:
        ref='t2t|hg38'
    shell:
        'python /gpfs/project/projects/medbioinf/projects/1000g-ont/extended_graph/svtig_to_single_contig.py -g {input.ref} -v {input.vcf}'

rule rename_pav_output:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav_{ref}/pav_svtigs.vcf.gz_merged.vcf'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav_{ref}/pav_svtigs_merged.vcf'
    wildcard_constraints:
        ref='t2t|hg38'
    resources:
        mem_total_mb=5000
    shell:
        'mv {input} {output}'