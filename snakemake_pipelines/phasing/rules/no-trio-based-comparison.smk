configfile: '../config.yaml'
flag_indel={'snp': '', 'indel': '--indels'}
wh_compare_flag={'snp': '--only-snvs', 'indel': ''}

rule no_trio_comparision:
    input:
        longread='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/longread/{sample}/{vtype}.vcf',
        nygc=lambda wildcards: '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-phased/no-trios/{sample}/filtered.vcf'
    output:
        pairwise='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/no-trio-comparision/{sample}/{vtype}/pairwise.tsv',
        multiway='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/no-trio-comparision/{sample}/{vtype}/multiway.tsv'
    wildcard_constraints:
        sample='(?:NA|HG)\d{5}',
        vtype='[a-z]{3,5}'
    conda:
        'envs/whatshap.yaml'
    params:
        indel_flag=lambda wildcards: wh_compare_flag[wildcards.vtype]
    resources:
        runtime_hrs=lambda wildcards, attempt: 3 * attempt,
        mem_total_mb=lambda wildcards, attempt: 10240 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 10240 * attempt
    shell:
        '''
        whatshap compare --sample {wildcards.sample} {params.indel_flag} --names longread,nygc --tsv-pairwise {output.pairwise} --tsv-multiway {output.multiway} {input.longread} {input.nygc}
        '''
