configfile: 'config.yaml'
flag_indel={'snp': '', 'indel': '--indels'}
wh_compare_flag={'snp': '--only-snvs', 'indel': ''}

rule extract_families:
    input:
        vcf='/gpfs/project/projects/medbioinf/users/spani/files/vcf/1000GP/NYGC/genotyped/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chr}.recalibrated_variants.vcf.gz'
    params:
        samples=lambda wildcards: ','.join((wildcards.family).split('_'))
    output:
        vcf=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-genotypes/trios/{family}/{chr}.vcf.gz'),
        vcf_index=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-genotypes/trios/{family}/{chr}.vcf.gz.tbi')
    conda:
        'envs/preprocessing.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: 5 * attempt,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 512 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 512 * attempt
    shell:
        '''
        bcftools view -s {params.samples} {input.vcf} | bgzip > {output.vcf}
        bcftools index --tbi -o {output.vcf_index} {output.vcf}
        '''

rule filter_families:
    input:
        vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-genotypes/trios/{family}/{chr}.vcf.gz',
        vcf_index='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-genotypes/trios/{family}/{chr}.vcf.gz.tbi'
    output:
        vcf=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-genotypes/trios/{family}/{chr}_filtered.vcf')
    conda:
        'envs/preprocessing.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: attempt,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 512 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 512 * attempt
    shell:
        '''
        bcftools view -f .,PASS --trim-alt-alleles -c 1 {input.vcf} | bcftools annotate -x INFO,FORMAT > {output.vcf}
        '''

rule extract_samples:
    input:
        vcf='/gpfs/project/projects/medbioinf/users/spani/files/vcf/1000GP/NYGC/genotyped/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chr}.recalibrated_variants.vcf.gz'
    output:
        vcf=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-genotypes/no-trios/{sample}/{chr}.vcf.gz'),
        vcf_index=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-genotypes/no-trios/{sample}/{chr}.vcf.gz.tbi')
    conda:
        'envs/preprocessing.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: 5 * attempt,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 512 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 512 * attempt
    shell:
        '''
        bcftools view -s {wildcard.sample} {input.vcf} | bgzip > {output.vcf}
        bcftools index --tbi -o {output.vcf_index} {output.vcf}
        '''

rule filter_samples:
    input:
        vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-genotypes/no-trios/{sample}/{chr}.vcf.gz',
        vcf_index=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-genotypes/no-trios/{sample}/{chr}.vcf.gz.tbi'),
    output:
        vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-genotypes/no-trios/{sample}/{chr}_filtered.vcf'
    conda:
        'envs/preprocessing.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: attempt,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 512 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 512 * attempt
    shell:
        '''
        bcftools view -f .,PASS --trim-alt-alleles -c 1 {input.vcf} | bcftools annotate -x INFO,FORMAT > {output.vcf}
        '''

rule extract_phased_trio:
    input:
        stat='/gpfs/project/projects/medbioinf/users/spani/files/vcf/1000GP/NYGC/phased/1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'
    params:
        samples=lambda wildcards: ','.join((wildcards.family).split('_'))
    output:
        stat=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-phased/trios/{family}/{chr}.vcf.gz')
    conda:
        'envs/preprocessing.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: 3 * attempt,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 512 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 512 * attempt
    shell:
        '''
        bcftools view -s {params.samples} -f .,PASS --trim-alt-alleles -c 1 {input.stat} | bcftools annotate -x INFO,FORMAT > {output.stat}
        '''

rule concat_stat_phased_trio:
    input:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-phased/trios/{{family}}/{chr}_filtered.vcf', chr=config['chromosome'])
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-phased/trios/{family}/filtered.vcf'
    conda:
        'envs/preprocessing.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: attempt,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        '''
        bcftools concat -o {output} {input}
        '''

rule extract_phased_single:
    input:
        stat='/gpfs/project/projects/medbioinf/users/spani/files/vcf/1000GP/NYGC/phased/1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'
    output:
        stat=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-phased/no-trios/{sample}/{chr}.vcf.gz')
    conda:
        'envs/preprocessing.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: 3 * attempt,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 512 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 512 * attempt
    shell:
        '''
        bcftools view -s {wildcard.sample} -f .,PASS --trim-alt-alleles -c 1 {input.stat} | bcftools annotate -x INFO,FORMAT > {output.stat}
        '''

rule concat_stat_phased_single:
    input:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-phased/no-trios/{{sample}}/{chr}_filtered.vcf', chr=config['chromosome'])
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-phased/no-trios/{sample}/filtered.vcf'
    conda:
        'envs/preprocessing.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: attempt,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        '''
        bcftools concat -o {output} {input}
        '''