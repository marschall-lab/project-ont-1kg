flag_indel={'snp': '', 'indel': '--indels'}
wh_compare_flag={'snp': '--only-snvs', 'indel': ''}

rule phase_trio:
    input:
        vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-genotypes/trios/{family}/{chr}_filtered.vcf',
        ref='/gpfs/project/projects/medbioinf/users/spani/files/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta',
        ped='/gpfs/project/projects/medbioinf/users/spani/files/other/1000GP/pedigree.ped'
    output:
        vcf=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/trio_phase/{family}/{chr}_{vtype}.vcf'),
        out='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/trio_phase/{family}/{chr}_{vtype}.out'
    conda:
        '../envs/whatshap.yaml'
    params:
        indel_flag=lambda wildcards: flag_indel[wildcards.vtype]
    resources:
        runtime_hrs=lambda wildcards, attempt: 2,
        mem_total_mb=lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        '''
        whatshap phase -o {output.vcf} --chromosome {wildcards.chr} {params.indel_flag} -r {input.ref} --ped {input.ped} {input.vcf} 2>&1 | tee {output.out} 
        '''

rule concat_trio_phase:
    input:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/trio_phase/{{family}}/{chr}_{{vtype}}.vcf', chr=config['chromosome'])
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/trio_phase/{family}/{vtype}.vcf'
    wildcard_constraints:
        chr='[c][h][r][0-9X]{1,2}',
        sample='(?:NA|HG)\d{5}',
        vtype='[a-z]{3,5}'
    conda:
        '../envs/preprocessing.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: 3 * attempt,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 512 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 512 * attempt
    shell:
        '''
        bcftools concat -o {output} {input}
        '''

rule stats_trio_phase:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/trio_phase/{family}/{vtype}.vcf'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/trio_phase/{family}/{vtype}.stats.tsv'
    conda:
        '../envs/whatshap.yaml'
    resources:
        runtime_hrs=0,
        runtime_min=10
    shell:
        'whatshap stats --tsv={output} {input}'    


rule longread_trio_phase:
    input:
        vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-genotypes/trios/{family}/{chr}_filtered.vcf',
        ref='/gpfs/project/projects/medbioinf/users/spani/files/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta',
        bam= lambda wildcards: expand('/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/hg38/{sample}.hg38.cram', sample=wildcards.family.split('_')),
        ped='/gpfs/project/projects/medbioinf/users/spani/files/other/1000GP/pedigree.ped'
    output:
        vcf=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/longread_trio_phase/{family}/{chr}_{vtype}.vcf'),
        out='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/longread_trio_phase/{family}/{chr}_{vtype}.out'
    conda:
        '../envs/whatshap.yaml'
    params:
        indel_flag=lambda wildcards: flag_indel[wildcards.vtype]
    resources:
        runtime_hrs=lambda wildcards, attempt: attempt,
        mem_total_mb=lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        '''
        whatshap phase -o {output.vcf} --chromosome {wildcards.chr} {params.indel_flag} -r {input.ref} --ped {input.ped} {input.vcf} {input.bam} 2>&1 | tee {output.out} 
        '''

rule concat_longread_trio:
    input:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/longread_trio_phase/{{family}}/{chr}_{{vtype}}.vcf', chr=config['chromosome'])
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/longread_trio_phase/{family}/{vtype}.vcf'
    wildcard_constraints:
        chr='[c][h][r][0-9X]{1,2}',
        sample='(?:NA|HG)\d{5}',
        vtype='[a-z]{3,5}'
    conda:
        '../envs/preprocessing.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: 3 * attempt,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 512 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 512 * attempt
    shell:
        '''
        bcftools concat -o {output} {input}
        '''

rule stats_longread_trio:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/longread_trio_phase/{family}/{vtype}.vcf'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/longread_trio_phase/{family}/{vtype}.stats.tsv'
    conda:
        '../envs/whatshap.yaml'
    resources:
        runtime_hrs=0,
        runtime_min=10
    shell:
        'whatshap stats --tsv={output} {input}' 