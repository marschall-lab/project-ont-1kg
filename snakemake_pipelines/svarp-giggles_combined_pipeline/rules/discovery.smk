# run svarp
rule svarp:
    input:
        reads = '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz',
        haplotag = '/gpfs/project/projects/medbioinf/data/share/globus/1000g-ont/read_tags/{sample}/{sample}_read_tags.tsv',
        alignment = '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/gaf/{sample}.sorted.gaf.gz',
        gfa = '/gpfs/project/projects/medbioinf/users/spani/files/gfa/HengLi/chm13-90c.r518_tagged.gfa',
        vcf = '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/annotating-paths/vcf/chm13-90c.r518_filtered.vcf.gz'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_H1.fa',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_H2.fa',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_untagged.fasta'
    params:
        outdir = '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}'
    log:
        "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}.log"
    conda:
        "../envs/svarp.yml"
    resources:
        mem_total_mb=40000,
        runtime_hrs=24,
        runtime_min=1
    priority: 1
    shell:
        """
        /gpfs/project/projects/medbioinf/projects/1000g-ont/svarp/build/svarp -a {input.alignment} -g {input.gfa} --fasta {input.reads} -i {wildcards.sample} -d 500 -s 5 --phase {input.haplotag} --e -o {params.outdir} 2>&1 | tee {log}
        """