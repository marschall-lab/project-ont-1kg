#compress vcf
rule compress_vcf:
    input:
        "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{filename}.vcf"
    output:
        vcf="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{filename}.vcf.gz",
        tbi="/gpfs/project/projects/medbioinf/users/spani/results/1000GP/augmented_graph/{filename}.vcf.gz.tbi"
    shell:
        """
        bgzip -c {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """