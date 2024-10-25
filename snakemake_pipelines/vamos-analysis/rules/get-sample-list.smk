import subprocess
import re
import gzip
import pandas as pd

#Finding sample list used for 1000GP project by us
path=config['path_to_t2t_cram']
ls_command = 'ls '+path
process = subprocess.Popen(ls_command.split(), stdout=subprocess.PIPE)
ls_out, ls_err = process.communicate()
sample_re = re.compile('(?:NA|HG)\d{5}')

cram_sample_list = []
for s in ls_out.decode('utf-8').split('\n'):
    if not s.endswith('cram'):
        continue
    cram_sample_list.append(s.split('.')[0])

# NYGC list
df = pd.read_csv('resources/igsr_sample_data.tsv', sep="\t")
nygc_sample_list = df["Sample name"].to_list()

t2t_samples = sorted(list(set(cram_sample_list).intersection(nygc_sample_list)))