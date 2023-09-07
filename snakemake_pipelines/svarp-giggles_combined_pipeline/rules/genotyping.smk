
# run genotyping
rule giggles:
    input:
        reads='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz',
        reads_gzi='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz.gzi',
        reads_fai='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz.fai',
        haplotag='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/haplotagging-results/GRCh38/{sample}/{sample}.tsv',
        alignment='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/gaf/{sample}.sorted.gaf',
        alignment_index='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/gaf/{sample}.sorted.gaf.gai',
        gfa='/gpfs/project/projects/medbioinf/users/spani/files/gfa/HengLi/chm13-90c.r518_tagged.gfa',
        vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/vcf/panel.vcf.gz'
    output:
        vcf=temp("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/{sample}.vcf")
    log:
        stderr="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/{sample}.stderr",
        stdout="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/genotypes/{sample}.stdout"
    resources:
        mem_total_mb=50000,
        runtime_hrs=48,
        runtime_min=1
    priority: 1
    shell:
        """
        set +u
        source ~/.bashrc
        conda activate giggles-dev
        set -u
        giggles genotype --read-fasta {input.reads} --sample {wildcards.sample} --realign-mode edit -o {output.vcf} --rgfa {input.gfa} --haplotag-tsv {input.haplotag} {input.vcf} {input.alignment} > {log.stdout} 2> {log.stderr}
        """