import argparse
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import seaborn
import statistics
import sys

def is_close_to_reference(info_fields, pos):
    start = int(pos[1])
    end = int(info_fields['END'])
    ref_len = end-start
    rus = info_fields['RU'].split(',')
    vntr_rus = info_fields['ALTANNO_H1'].split(',')
    vntr_len = sum([len(rus[int(i)]) for i in vntr_rus])
    if ref_len == 0:
        return False
    if abs(ref_len-vntr_len)/ref_len < 0.1:
        return True
    else:
        return False

def parse_info(line):
    info_fields={}
    for field in line.split(';')[:-1]:
        assert '=' in field
        key, value = field.split('=')
        info_fields[key] = value
    return info_fields

def determine_pairs(h1, h2, o1, o2, mean):
    if (abs(h1-o1)+abs(h2-o2)) < (abs(h1-o2)+abs(h2-o1)):
        return (h1-mean, o1-mean), (h2-mean, o2-mean)
    else:
        return (h1-mean, o2-mean), (h2-mean, o1-mean)

def run(hgsvc, ont, sample, output):
    
    # reading hgsvc vntr calls
    hap1, hap2 = hgsvc.split(',')
    hap1_reader = open(hap1, 'r')
    hap2_reader = open(hap2, 'r')
    hgsvc_data = {}
    count = 0
    # reading haplotype 1
    while True:
        hap1_line = hap1_reader.readline()
        if not hap1_line:
            break
        if hap1_line[0] == '#':
            continue
        hap1_line = hap1_line.rstrip().split('\t')
        hgsvc_data[(hap1_line[0], hap1_line[1])] = [parse_info(hap1_line[7]), None]
    # reading haplotype 2
    while True:
        hap2_line = hap2_reader.readline()
        if not hap2_line:
            break
        if hap2_line[0] == '#':
            continue
        hap2_line = hap2_line.rstrip().split('\t')
        try:
            hgsvc_data[(hap2_line[0], hap2_line[1])][1] = parse_info(hap2_line[7])
            count += 1
        except KeyError:
            continue
    hap1_reader.close()
    hap2_reader.close()

    print(f'Number of positions where HGSVC has VNTRs for both haplotypes: {count}', file=sys.stderr)

    # reading ONT vntr calls
    ont_reader = open(ont, 'r')
    ont_data = {}
    pos_list = []
    count = 0
    while True:
        line = ont_reader.readline()
        if not line:
            break
        if line[0] == '#':
            continue
        count += 1
        line = line.rstrip().split('\t')
        if (line[0], line[1]) in hgsvc_data:
            if hgsvc_data[(line[0], line[1])][1] == None:
                continue
            ont_data[(line[0], line[1])] = parse_info(line[7])
            pos_list.append((line[0], line[1]))
    ont_reader.close()

    print(f'Number of positions where ONT has VNTRs for both haplotypes: {count}', file=sys.stderr)

    # plotting the density of RU counts of VNTRs
    hgsvc_counts = []
    ont_counts = []
    hgsvc_counts_subset100 = []
    ont_counts_subset100 = []
    
    for pos in pos_list:
        hgsvc_count_1 = int(hgsvc_data[pos][0]['LEN_H1'])
        is_hg1_ref = is_close_to_reference(hgsvc_data[pos][0], pos)
        hgsvc_count_2 = int(hgsvc_data[pos][1]['LEN_H1'])
        is_hg2_ref = is_close_to_reference(hgsvc_data[pos][1], pos)
        ont_count_1 = int(ont_data[pos]['LEN_H1'])
        try:
            ont_count_2 = int(ont_data[pos]['LEN_H2'])
        except KeyError:
            ont_count_2 = int(ont_data[pos]['LEN_H1'])
        ru_lens = [len(i) for i in ont_data[pos]['RU'].split(',')]
        ru_len_mean = statistics.fmean(ru_lens)
        # mean number of repeat units in the reference
        n_rus = (int(ont_data[pos]['END'])-int(pos[1]))/ru_len_mean
        
        pair1, pair2 = determine_pairs(hgsvc_count_1, hgsvc_count_2, ont_count_1, ont_count_2, n_rus)
        #if pair1[0] > 100 or pair1[1] > 100 or pair2[0] > 100 or pair2[1] > 100:
        #    continue
        # only considering pairs where the hgsvc VNTRs have different RU numbers than the ref (threshold of 0.5)
        #if abs(pair1[0]) > 1:
        if not is_hg1_ref:
            hgsvc_counts.append(pair1[0])
            ont_counts.append(pair1[1])
            if abs(pair1[0]) < 101 and abs(pair1[1]) < 101:
                hgsvc_counts_subset100.append(pair1[0])
                ont_counts_subset100.append(pair1[1])
        #if abs(pair2[0]) > 1:
        if not is_hg2_ref:
            hgsvc_counts.append(pair2[0])
            ont_counts.append(pair2[1])
            if abs(pair2[0]) < 101 and abs(pair2[1]) < 101:
                hgsvc_counts_subset100.append(pair2[0])
                ont_counts_subset100.append(pair2[1])

    if '_' in sample:
        hgsvc_sample, ont_sample = sample.split('_')
    else:
        hgsvc_sample = sample
        ont_sample = sample
    
    stats_writer = open(f'{output}/{sample}.stats', 'w')

    print(f'Number of positions plotted: {len(hgsvc_counts)}', file=sys.stderr)
    g = seaborn.JointGrid(x=hgsvc_counts, y=ont_counts, height=10, ratio=5)
    g.plot_marginals(seaborn.histplot, bins=100)
    g.plot_joint(plt.hexbin, bins='log', gridsize=50, cmap='Blues')
    g.set_axis_labels("Difference of RU Counts in HGSVC VNTRs v/s Reference", "Difference of RU Counts in ONT VNTRs v/s Reference", fontsize=20)
    g.fig.suptitle(f'HGSVC: {hgsvc_sample}    ONT: {ont_sample}\nPCC: {pearsonr(hgsvc_counts, ont_counts)[0]}', fontsize=20)
    cbar_ax = g.fig.add_axes([.95, .2, .02, .6])  # x, y, width, height
    plt.colorbar(cax=cbar_ax)
    g.savefig(f'{output}/{sample}-RU.svg')
    print(f'01\tPCC_RU-Scatter_Whole\t{pearsonr(hgsvc_counts, ont_counts)[0]}', file=stats_writer)

    f = seaborn.JointGrid(x=hgsvc_counts_subset100, y=ont_counts_subset100, height=10, ratio=5, ylim=(-100, 100), xlim=(-100,100))
    #g = seaborn.JointGrid(x=hgsvc_counts, y=ont_counts, height=10, ratio=5)
    f.plot_marginals(seaborn.histplot, bins=range(0, 100))
    f.plot_joint(plt.hexbin, bins='log', gridsize=50, cmap='Blues')
    f.set_axis_labels("Difference of RU Counts in HGSVC VNTRs v/s Reference", "Difference of RU Counts in ONT VNTRs v/s Reference", fontsize=20)
    f.fig.suptitle(f'HGSVC: {hgsvc_sample}    ONT: {ont_sample}\nPCC: {pearsonr(hgsvc_counts_subset100, ont_counts_subset100)[0]}', fontsize=20)
    cbar_ax = f.fig.add_axes([.95, .2, .02, .6])  # x, y, width, height
    plt.colorbar(cax=cbar_ax)
    f.savefig(f'{output}/subset100-{sample}-RU.svg')
    print(f'02\tPCC_RU-Scatter_Subset100\t{pearsonr(hgsvc_counts_subset100, ont_counts_subset100)[0]}', file=stats_writer)

    stats_writer.close()
    
    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='plot-hgsvc-qc.py', description="Creates plots for HGSVC-based QC")
    parser.add_argument("-hgsvc", required=True, help="Comma-separated <hap1,hap2> list of VNTRS from HGSVC")
    parser.add_argument("-ont", required=True, help="ONT Vamos VCF")
    parser.add_argument("-sample", required=True, help="Sample name")
    parser.add_argument("-output", required=True, help="Output prefix")

    options = parser.parse_args()

    run(**vars(options))