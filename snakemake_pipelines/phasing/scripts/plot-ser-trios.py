import argparse
import pandas
import numpy
import matplotlib.pyplot as plt
import matplotlib
from collections import namedtuple

def read_tsv(file):
    
    tsv = pandas.read_csv(file, sep='\t', header=0)
    sample = tsv['#sample'].to_numpy()[0]
    tsv = tsv[['chromosome', 'dataset_name0', 'dataset_name1', 'all_assessed_pairs', 'all_switches', 'blockwise_hamming']]
    data = {
        ('trio', 'longread'): numpy.array([0, 0, 0]),
        ('trio', 'trio-longread'): numpy.array([0, 0, 0]),
        ('trio', 'nygc'): numpy.array([0, 0, 0]),
        ('longread', 'trio-longread'): numpy.array([0, 0, 0]),
        ('longread', 'nygc'): numpy.array([0, 0, 0]),
        ('trio-longread', 'nygc'): numpy.array([0, 0, 0]),
        }
    
    for line in tsv.to_numpy():
        _, d1, d2, x, y, z = line
        data[(d1, d2)] += numpy.array([x, y, z])
    
    return sample, data

def read_metadata(file):
    metadata = pandas.read_csv(file, sep='\t', header=0)
    metadata = metadata[["Sample name", "Population code", "Superpopulation code"]]
    sample2pop = {}
    pop2color = {'AFR': '#FFCD33', 'AMR': '#FF3D3D', 'EAS': '#ADFF3D', 'EUR': '#64EBFF', 'SAS': '#FF30FF'}
    for line in metadata.to_numpy():
        sample, _, pop = line
        sample2pop[sample] = pop
    
    return sample2pop, pop2color

parser = argparse.ArgumentParser(prog='plot-ser-boxplot.py', description="Plotting the phasing comparison results")
parser.add_argument('-tsvs', metavar='TSVS', help='File containing list of Pairwise TSV output files from whatshap compare')
parser.add_argument('-output', metavar='OUTPUT', help='output directory')
args = parser.parse_args()

families = {
  '2418': 'NA19818_NA19819_NA19828'.split('_'),
  'CLM16': 'HG01256_HG01257_HG01258'.split('_'),
  'SH006': 'HG00418_HG00419_HG00420'.split('_'),
  'Y077': 'NA19128_NA19127_NA19129'.split('_'),
}

sample2family = {}

for family, samples in families.items():
    for sample in samples:
        sample2family[sample] = family

tsv_list = pandas.read_csv(args.tsvs, header=None).to_numpy()[:,0]

sample_results = {}
for tsv in tsv_list:
    sample, data = read_tsv(tsv)
    sample_results[sample] = data

# switch error rate against the nygc phasing
longread_ser = {}
trio_ser = {}
triolongread_ser = {}

for sample in sample2family.keys():
    data = sample_results[sample]
    longread_ser[sample] = data[('longread', 'nygc')][1]/data[('longread', 'nygc')][0]
    trio_ser[sample] = data[('trio', 'nygc')][1]/data[('trio', 'nygc')][0]
    triolongread_ser[sample] = data[('trio-longread', 'nygc')][1]/data[('trio-longread', 'nygc')][0]


## Plotting

fig = plt.figure(figsize =(10, 5))
x_pos = [1,2,3,5,6,7,9,10,11,13,14,15]
x_labels = list(sample2family.keys())

plt.plot(x_pos, [y for y in longread_ser.values()], 'bo', label='Long-read phasing')
plt.plot(x_pos, [y for y in triolongread_ser.values()], 'go', label = 'Long-read + Trio phasing')
plt.plot(x_pos, [y for y in trio_ser.values()], 'ro', label = 'Trio phasing')
for x in [4,8,12]:
    plt.axvline(x=x, color='black', ls='--', lw=1)
plt.ylim((0,0.020))
plt.xticks(x_pos, labels=x_labels, rotation=45)
plt.ylabel('Switch Error Rate', fontsize=15)
plt.title('WhatsHap compare against NYGC Statistical Phasing', fontsize=15)
plt.legend()
plt.tight_layout()
plt.savefig(args.output+'/trio-ser.svg')
plt.savefig(args.output+'/trio-ser.pdf')