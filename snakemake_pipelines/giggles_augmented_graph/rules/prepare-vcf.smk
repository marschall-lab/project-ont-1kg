configfile: 'config/config.yaml'

# creating the list of pseudohaplotypes available for all callsets
pseudohaplotypes = {}
for callset in config['callsets'].keys():
    pseudohaplotypes[callset] = {}
    if config['callsets'][callset]['pseudohaplotypes'] == None:
        break
    for line in open(config['callsets'][callset]['pseudohaplotypes'], 'r'):
        line = line.strip()
        name = line.split('/')[-1].split('.')[0]
        pseudohaplotypes[callset][name] = line

# extract HPRC assemblies from AGC
# TODO: Can add a rule for doing this extraction
# if assemblies have been deleted, there is a script in /gpfs/project/projects/medbioinf/users/spani/files/fasta/HPRC/Assemblies_Yr1/ to extract them.

# align HPRC assemblies
rule align_hprc_assemblies:
    input:
        ref='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa',
        assembly='/gpfs/project/projects/medbioinf/users/spani/files/fasta/HPRC/Assemblies_Yr1/{sample}.{haplotype}.fa'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/hprc_assembly_mappings/{sample}.{haplotype}.gaf'
    log:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/hprc_assembly_mappings/{sample}.{haplotype}.log'
    wildcard_constraints:
        sample='|'.join(config['hprc_samples'])
    resources:
        runtime_hrs=12,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 80000 * attempt
    threads: 8
    shell:
        '/gpfs/project/projects/medbioinf/users/spani/packages/minigraph/minigraph --vc -cx lr {input.ref} {input.assembly} -t {threads} > {output} 2> {log}'


# align pseudohaplotypes
rule align_pseudohaplotypes:
    input:
        ref='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa',
        assembly=lambda wildcards: pseudohaplotypes[wildcards.callset][wildcards.name]
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/pseudo_assembly_mappings/{name}.gaf'
    log:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/pseudo_assembly_mappings/{name}.log'
    resources:
        runtime_hrs=12,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 80000 * attempt
    threads: 8
    shell:
        '/gpfs/project/projects/medbioinf/users/spani/packages/minigraph/minigraph --vc -cx lr {input.ref} {input.assembly} -t {threads} > {output} 2> {log}'


# make a list of HPRC assembly GAFs
rule list_hprc_assembly_mappings:
    input:
        expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{{callset}}/hprc_assembly_mappings/{sample}.{haplotype}.gaf', sample=config['hprc_samples'], haplotype=['1', '2'])
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/hprc_assembly_mappings/pathlist.txt'
    run:
        f = open(output[0], 'w')
        for name in input:
            print(name, file=f)
        f.close()  


# make a list of pseudo-haplotype GAFs
rule list_pseudohaplotype_assembly_mappings:
    input:
        lambda wildcards: expand('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{{callset}}/pseudo_assembly_mappings/{name}.gaf', name=[p for p in pseudohaplotypes[wildcards.callset].keys()])
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/pseudo_assembly_mappings/pathlist.txt'
    run:
        f = open(output[0], 'w')
        for name in input:
            print(name, file=f)
        f.close()  


# script to make the VCF taking the two lists made above
rule assemblies_to_vcf:
    input:
        ref='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa',
        hprc_path='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/hprc_assembly_mappings/pathlist.txt',
        pseudo_path='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/pseudo_assembly_mappings/pathlist.txt'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/vcf/{callset}.vcf'
    log:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/vcf/{callset}.log'
    conda:
        '../envs/basic.yml'
    resources:
        runtime_hrs=1,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 10*1024 * attempt
    shell:
        'python scripts/assembly-to-vcf.py -gfa {input.ref} -hprc-list {input.hprc_path} -pseudo-list {input.pseudo_path} --output {output} 2> {log}'


# filter the VCF
rule filter_vcf:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/vcf/{callset}.vcf'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/giggles_augmented_graph/{callset}/vcf/{callset}_filtered.vcf.gz'
    conda:
        '../envs/basic.yml'
    resources:
        runtime_hrs=1,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 2*1024 * attempt
    shell:
        '''
        bcftools view -f "PASS" -Oz -o {output} {input}
        bcftools index {output}
        '''