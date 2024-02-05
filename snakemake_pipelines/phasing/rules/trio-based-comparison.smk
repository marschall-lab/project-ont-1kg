configfile: '../config.yaml'
flag_indel={'snp': '', 'indel': '--indels'}
wh_compare_flag={'snp': '--only-snvs', 'indel': ''}

sample2family={}
for family in config['families']:
    for sample in family.split('_'):
        sample2family[sample] = family

rule trio_comparision:
    input:
        trio=lambda wildcards: '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/trio_phase/'+sample2family[wildcards.sample]+'/{vtype}.vcf',
        longread_trio=lambda wildcards: '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/longread_trio_phase/'+sample2family[wildcards.sample]+'/{vtype}.vcf',
        longread='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/phased-vcf/longread/{sample}/{vtype}.vcf',
        nygc=lambda wildcards: '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/data/nygc-phased/trios/'+sample2family[wildcards.sample]+'/filtered.vcf'
    output:
        comp12='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/trio-comparision/{sample}/{vtype}/trio_longread.txt',
        comp13='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/trio-comparision/{sample}/{vtype}/trio_triolongread.txt',
        comp1stat='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/trio-comparision/{sample}/{vtype}/trio_stat.txt',
        pairwise='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/trio-comparision/{sample}/{vtype}/pairwise.tsv',
        multiway='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/trio-comparision/{sample}/{vtype}/multiway.tsv'
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
        whatshap compare --sample {wildcards.sample} {params.indel_flag} --names trio,longread {input.trio} {input.longread} > {output.comp12}
        whatshap compare --sample {wildcards.sample} {params.indel_flag} --names trio,trio-longread {input.trio} {input.longread_trio} > {output.comp13}
        whatshap compare --sample {wildcards.sample} {params.indel_flag} --names trio,nygc {input.trio} {input.nygc} > {output.comp1stat}
        whatshap compare --sample {wildcards.sample} {params.indel_flag} --names trio,longread,trio-longread,nygc --tsv-pairwise {output.pairwise} --tsv-multiway {output.multiway} {input.trio} {input.longread} {input.longread_trio} {input.nygc}
        '''
