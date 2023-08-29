import subprocess
import re
import gzip
import pandas as pd

#Finding sample list used for 1000GP project by us
path='/gpfs/project/projects/medbioinf/data/share/globus/1000g-ont/hg38/'
ls_command = 'ls '+path
process = subprocess.Popen(ls_command.split(), stdout=subprocess.PIPE)
ls_out, ls_err = process.communicate()
sample_re = re.compile('(?:NA|HG)\d{5}')

cram_sample_list = []
for s in ls_out.decode('utf-8').split('\n'):
    if not s.endswith('cram'):
        continue
    cram_sample_list.append(s.split('.')[0])
    
vcf_sample_list = None
with gzip.open("/gpfs/project/projects/medbioinf/users/spani/results/1000GP/haplotagging-results/nygc-merged.vcf.gz", "rt") as f:
    for line in f:
        if line[0:2] == "##":
            continue
        vcf_sample_list = line.rstrip().split("\t")[9:]
        break

samples = list(set(cram_sample_list).intersection(vcf_sample_list))