import subprocess

#Finding sample list used for 1000GP project by us
path=config['path_to_ss_delly']
ls_command = 'ls '+path
process = subprocess.Popen(ls_command.split(), stdout=subprocess.PIPE)
ls_out, ls_err = process.communicate()

sample_list = []
for s in ls_out.decode('utf-8').split('\n'):
    if not s.endswith('bcf'):
        continue
    sample_list.append(s.split('.')[0])

rule convert_bcf_to_vcf:
    input:
        config['path_to_ss_delly']+'{sample}.bcf'
    output:
        'results/ss_delly_vcfs/{sample}.vcf'
    conda:
        '../envs/analysis.yml'
    shell:
        'bcftools view {input} > {output}'

rule convert_all:
    input:
        expand('results/ss_delly_vcfs/{sample}.vcf', sample=sample_list)
    output:
        'results/ss_delly_vcfs/done.chk'
    shell:
        'echo "All done" > {output}'