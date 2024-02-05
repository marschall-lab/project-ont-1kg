configfile: './config.yaml'
flag_indel={'snp': '', 'indel': '--indels'}
wh_compare_flag={'snp': '--only-snvs', 'indel': ''}

rule phase_sample_longread:
    input:
        vcf='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-genotypes/no-trios/{sample}/{chr}_filtered.vcf',
        bam='/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/cram/{sample}/alignments/{sample}.cram',
        ref='/gpfs/project/projects/medbioinf/users/spani/files/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta'
    output:
        vcf=temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/longread/{sample}/{chr}_{vtype}.vcf'),
        out='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/longread/{sample}/{chr}_{vtype}.out'
    params:
        indel_flag=lambda wildcards: flag_indel[wildcards.vtype]
    conda:
        'envs/whatshap.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: attempt,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        '''
        whatshap phase -o {output.vcf} --chromosome {wildcards.chr} --sample {wildcards.sample} {params.indel_flag} -r {input.ref} {input.vcf} {input.bam} 2>&1 | tee {output.out} 
        '''

rule concat_longread:
    input:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/longread/{{sample}}/{chr}_{{vtype}}.vcf', chr=config['chromosome'])
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/longread/{sample}/{vtype}.vcf'
    wildcard_constraints:
        chr='[c][h][r][0-9X]{1,2}',
        sample='(?:NA|HG)\d{5}',
        vtype='[a-z]{3,5}'
    conda:
        'envs/preprocessing.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: 3 * attempt,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 512 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 512 * attempt
    shell:
        '''
        bcftools concat -o {output} {input}
        '''