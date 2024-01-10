import subprocess
import re
import gzip
import pandas as pd

#Finding sample list used for 1000GP project by us
path='/gpfs/project/projects/medbioinf/data/share/globus/hhu-1000g-ont/hg38/'
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

chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']