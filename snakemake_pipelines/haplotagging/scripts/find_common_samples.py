import subprocess
import re
import gzip

#Finding sample list used for 1000GP project by us
path='/gpfs/project/projects/medbioinf/data/share/globus/1000g-ont/GRCh38/'
ls_command = 'ls '+path
process = subprocess.Popen(ls_command.split(), stdout=subprocess.PIPE)
ls_out, ls_err = process.communicate()
sample_re = re.compile('(?:NA|HG)\d{5}')

cram_sample_list = [s for s in ls_out.decode('utf-8').split('\n') if sample_re.match(s)]

vcf_sample_list = None
with gzip.open("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/haplotagging-results/nygc-merged.vcf.gz", "rt") as f:
    for line in f:
        if line[0:2] == "##":
            continue
        vcf_sample_list = line.rstrip().split("\t")[9:]
        break

common_samples = list(set(cram_sample_list).difference(vcf_sample_list))
for i in cram_sample_list:
    print(i)
