from matplotlib_venn import venn3, venn3_unweighted, venn3_circles
import matplotlib.pyplot as plt

nygc = '/home/samarendra/Work/genotyping-results-localcopy/3202-nygc-sample-list.tsv'
nygc_unrelated = '/home/samarendra/Work/genotyping-results-localcopy/2504_unrelated-sample-list.tsv'
vienna_ont = '/home/samarendra/Work/genotyping-results-localcopy/1kg-vienna-sample-list.tsv'
hprc = '/home/samarendra/Work/genotyping-results-localcopy/47-hprc-sample-list.txt'
hgsvc = '/home/samarendra/Work/genotyping-results-localcopy/44-hgsvc-phase2-sample-list.tsv'

nygc_samples = set()
nygc_unrelated_samples = set()
vienna_ont_samples = set()
hprc_samples = set()
hgsvc_samples = set()

with open(nygc, 'r') as file:
    for line in file:
        nygc_samples.add(line.rstrip())

with open(nygc_unrelated, 'r') as file:
    for line in file:
        nygc_unrelated_samples.add(line.rstrip())

with open(vienna_ont, 'r') as file:
    for line in file:
        vienna_ont_samples.add(line.rstrip())

with open(hprc, 'r') as file:
    for line in file:
        hprc_samples.add(line.rstrip())

with open(hgsvc, 'r') as file:
    for line in file:
        hgsvc_samples.add(line.rstrip())

plt.figure(figsize=(10,10))
v = venn3([nygc_samples, nygc_unrelated_samples, vienna_ont_samples], ['NYGC Complete Sample List', 'NYGC Unrelated Sample List', '1000GP-ONT Vienna Sample List'])
c = venn3_circles([nygc_samples, nygc_unrelated_samples, vienna_ont_samples])
plt.savefig('/home/samarendra/Work/genotyping-results-localcopy/sample-venn/venn1.svg')
plt.close()

vienna_ont_graph_set = vienna_ont_samples & nygc_samples
plt.figure(figsize=(10,10))
v = venn3([vienna_ont_graph_set, hprc_samples, hgsvc_samples], ['Working Set of 1000GP-ONT Vienna Sample List', 'HPRC Sample List', 'HGSVC Phase 2 Sample List'])
c = venn3_circles([vienna_ont_graph_set, hprc_samples, hgsvc_samples])
plt.savefig('/home/samarendra/Work/genotyping-results-localcopy/sample-venn/venn2.svg')
plt.close()