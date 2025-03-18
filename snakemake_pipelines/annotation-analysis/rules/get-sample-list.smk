import subprocess
import re
import gzip
import pandas as pd

#Finding sample list used for 1000GP project by us
path=config['path_to_t2t_cram']
ls_command = 'ls '+path
process = subprocess.Popen(ls_command.split(), stdout=subprocess.PIPE)
ls_out, ls_err = process.communicate()

cram_sample_list = []
for s in ls_out.decode('utf-8').split('\n'):
    if not s.endswith('cram'):
        continue
    cram_sample_list.append(s.split('.')[0])

# NYGC list
df = pd.read_csv('resources/igsr_sample_data.tsv', sep="\t")
nygc_sample_list = df["Sample name"].to_list()

t2t_samples = sorted(list(set(cram_sample_list).intersection(nygc_sample_list)))

# Finding the samples common to HGSVC
path=config['path_to_hgsvc_alignments']
ls_command = 'ls '+path
process = subprocess.Popen(ls_command.split(), stdout=subprocess.PIPE)
ls_out, ls_err = process.communicate()

hgsvc_sample_list_all = []
for s in ls_out.decode('utf-8').split('\n'):
    if not s.endswith('.vrk-ps-sseq.asm-hap1.t2tv2.sort.bam'):
        continue
    hgsvc_sample_list_all.append(s.split('/')[-1][0:7])

# set of samples intersecting with our callset
hgsvc_samples = sorted(list(set(t2t_samples).intersection(hgsvc_sample_list_all)))

# all hgsvc3 assembly samples except HG00514 (asked to be removed by P. Ebert)
hgsvc_sample_list_all.remove('HG00514')