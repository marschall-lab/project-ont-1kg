
# run svarp
rule svarp:
    input:
        reads='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz',
        reads_gzi='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz.gzi',
        reads_fai='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz.fai',
        haplotag='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/haplotagging-results/GRCh38/{sample}/{sample}.tsv',
        alignment='/gpfs/project/projects/medbioinf/data/share/globus/1000g-ont/gaf/{sample}.gaf.gz',
        gfa='/gpfs/project/projects/medbioinf/users/spani/files/gfa/HengLi/chm13-90c.r518_tagged.gfa'
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
        mem_total_mb=196000,
        runtime_hrs=96,
        runtime_min=1
    priority: 1
    shell:
        '''
        module load Python/3.11.3
        module load Minimap2/2.17
        module load SamTools/1.6
        module load gcc/10.2.0
        export PATH=$PATH:/gpfs/project/projects/medbioinf/users/spani/packages/wtdbg2
        /gpfs/project/projects/medbioinf/projects/1000g-ont/svarp/build/svarp -a {input.alignment} -g {input.gfa} --fasta {input.reads} -i {wildcards.sample} -d 500 -s 5 --phase {input.haplotag} --e -o {params.outdir} 2> {log.stderr} 1> {log.stdout}
        module unload gcc/10.2.0
        module unload Minimap2/2.17
        module unload Python/3.11.3
        module unload SamTools/1.6
        '''