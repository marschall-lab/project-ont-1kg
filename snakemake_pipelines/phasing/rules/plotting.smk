rule list_trios:
    input:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/trio-comparision/{sample}/snp/pairwise.tsv', sample=family_sample_list)
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/trio-comparision/pairwise-list.tsv'
    resources:
        runtime_hrs=0,
        runtime_min=10
    run:
        f = open(output[0], 'w')
        for name in input:
            print(name, file=f)
        f.close()

rule plot_trios:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/trio-comparision/pairwise-list.tsv'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/plots/trio-ser.svg',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/plots/trio-ser.pdf'
    params:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/plots/'
    conda:
        '../envs/basic.yaml'
    shell:
        'python scripts/plot-ser-trios.py -tsvs {input} -output {params}'


rule list_non_trios:
    input:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/no-trio-comparision/{sample}/snp/pairwise.tsv', sample=samples)
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/no-trio-comparision/pairwise-list.tsv'
    resources:
        runtime_hrs=0,
        runtime_min=10
    run:
        f = open(output[0], 'w')
        for name in input:
            print(name, file=f)
        f.close()

rule plot_non_trios:
    input:
        tsv='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/no-trio-comparision/pairwise-list.tsv',
        meta='/gpfs/project/projects/medbioinf/users/spani/files/other/1000GP/igsr_sample_data.tsv'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/plots/no-trio-ser.svg',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/plots/no-trio-ser.pdf'
    params:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/phasing-results/plots/'
    conda:
        '../envs/basic.yaml'
    shell:
        'python scripts/plot-ser-boxplot.py -tsvs {input.tsv} -meta {input.meta} -output {params}'   

        