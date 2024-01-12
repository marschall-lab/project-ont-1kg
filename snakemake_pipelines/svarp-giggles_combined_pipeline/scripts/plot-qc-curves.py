import random
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerTuple
import pandas
from cyvcf2 import VCF
import argparse
import numpy as np
from collections import namedtuple, Counter

def process_SV(name):
    _, _, type, bub, len = name.split('-')
    return type, bub, len

def plot_audano_curve(svs_per_sample, metadata, output):
    # TODO: Create Legend. Figure out how to better represent the superpopulations when number of samples increases
    color_by_pop_1 = {'AFR': '#ffd845', 'AMR': '#bfbcb0', 'EAS': '#bfbcb0', 'EUR': '#bfbcb0', 'SAS': '#bfbcb0'}
    samples = list(svs_per_sample.keys())
    # create legend
    handles = []
    labels = []
    handles.append([mpatches.Patch(facecolor=c, label='Singletons (AC = 1)') for c in ['#ffd845', '#bfbcb0']])
    labels.append('Singletons (AC = 1)')
    handles.append([mpatches.Patch(facecolor=c, label='Doubletons (AC = 2)') for c in ['#ffca04', '#96917d']])
    labels.append('Doubletons (AC = 2)')
    handles.append([mpatches.Patch(facecolor=c, label='Polymorphic SVs (AC > 2)') for c in ['#c29a00', '#96917d']])
    labels.append('Polymorphic SVs (AC > 2)')
    handles.append([mpatches.Patch(facecolor=c, label='Major SVs (AF >= 50%)') for c in ['#826600', '#514e42']])
    labels.append('Major SVs (AF >= 50%)')
    handles.append(mpatches.Patch(facecolor='red', label='Shared SVs (AF = 100%)'))
    labels.append('Shared SVs (AF = 100%)')
    colors_1 = []
    colors_2 = []
    colors_3 = []
    colors_4 = []
    for s in samples:
        colors_1.append('#ffd845' if metadata[metadata['Sample name'] == s]['Superpopulation code'].values[0] == 'AFR' else '#bfbcb0')
        colors_2.append('#ffca04' if metadata[metadata['Sample name'] == s]['Superpopulation code'].values[0] == 'AFR' else '#9e9a87')
        colors_3.append('#c29a00' if metadata[metadata['Sample name'] == s]['Superpopulation code'].values[0] == 'AFR' else '#7a7563')
        colors_4.append('#826600' if metadata[metadata['Sample name'] == s]['Superpopulation code'].values[0] == 'AFR' else '#514e42')
    count_shared_svs = []
    count_major_svs = []
    count_polymorphic_svs = []
    count_doubletons_svs = []
    count_singleton_svs = []
    sv_counter = Counter()
    for n, sample in enumerate(samples):
        sample_svs = [sv.name for sv in svs_per_sample[sample] if int(sv.len) > 49]
        sv_counter.update(sample_svs)
        shared_counter = 0
        major_counter = 0
        polymorphic_counter = 0
        doubleton_counter = 0
        singleton_counter = 0
        for value in sv_counter.values():
            if value == n+1:
                singleton_counter += 1
                doubleton_counter += 1
                polymorphic_counter += 1
                major_counter += 1
                shared_counter += 1
            elif value >= (n+1)/2:
                singleton_counter += 1
                doubleton_counter += 1
                polymorphic_counter += 1
                major_counter += 1
            elif value > 2:
                singleton_counter += 1
                doubleton_counter += 1
                polymorphic_counter += 1
            elif value == 2:
                singleton_counter += 1
                doubleton_counter += 1
            elif value == 1:
                singleton_counter += 1
        count_shared_svs.append(shared_counter)
        count_major_svs.append(major_counter)
        count_polymorphic_svs.append(polymorphic_counter)
        count_doubletons_svs.append(doubleton_counter)
        count_singleton_svs.append(singleton_counter)
    
    fig, ax = plt.subplots(figsize = (20,10))
    ax.bar(samples, count_singleton_svs, color=colors_1)
    ax.bar(samples, count_doubletons_svs, color=colors_2)
    ax.bar(samples, count_polymorphic_svs, color=colors_3)
    ax.bar(samples, count_major_svs, color=colors_4)
    ax.bar(samples, count_shared_svs, color='red')
    ax.grid(axis='y', color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.ylabel("Cumulative SVs", fontsize=15)
    plt.xlabel("Samples", fontsize=15)
    plt.yticks(fontsize=10)
    plt.figlegend(handles, labels, framealpha=1, frameon=True, loc='upper left', bbox_to_anchor=(0.1, 0.99), handler_map = {list: HandlerTuple(None)})
    plt.tight_layout()
    plt.savefig('%s/audano-curve.png'%(output))
    plt.close()
    '''
    # Removing all the singletons and replotting.
    singleton_svs_list = []
    for sv_name, count in sv_counter.items():
        if count == 1:
            singleton_svs_list.append(sv_name)
    
    fig, ax = plt.subplots(figsize = (20,10))
    ax.bar(samples, count_cumulative_svs, label='cumulative svs', color=colors_1)
    ax.bar(samples, count_shared_svs, label='shared svs', color='black')
    ax.grid(axis='y', color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.ylabel("Cumulative SVs", fontsize=15)
    plt.xlabel("Samples", fontsize=15)
    plt.yticks(fontsize=10)
    plt.figlegend(handles, labels, framealpha=1, frameon=True, loc='upper left', bbox_to_anchor=(0.1, 0.9))
    plt.tight_layout()
    plt.savefig('%s/audano-curve_no-singletons.png'%(output))
    plt.close()
    '''
    
def plot_sv_size_distribution(size_dist, output):
    
    def construct_plot(s, r, o, b, offset):
        mean_list = []
        sum_list = []
        for sv_type in ['DEL', 'INS', 'COMPLEX']:
            data = s[sv_type]
            mean = []
            sum = []
            stddev = []
            for bin_data in data:
                tmp = np.array(bin_data)
                mean.append(np.mean(tmp))
                sum.append(np.sum(tmp))
                stddev.append(tmp.std())
            mean_list.append(mean)
            sum_list.append(sum)
            x_postion = np.array(range(len(data)))
            x_labels = ((x_postion+1)*b)+offset
            width = 0.8
            # Plotting average count over all samples
            fig = plt.figure(figsize =(10, 10))
            plt.bar(x_postion, mean, width=width, edgecolor='black')
            #plt.errorbar(x_postion, mean, yerr=stddev, fmt='none', ecolor='Black', elinewidth=2)
            plt.xlabel('Size (in bp)', fontsize=15)
            plt.ylabel('Average Count per Sample', fontsize=15)
            plt.xticks(x_postion+0.5, x_labels, rotation=60, fontsize=5)
            plt.yticks(fontsize=10)
            plt.title(sv_type, fontsize=20)
            plt.tight_layout()
            plt.savefig('%s/average_count-size-distribution_%s_%s.png'%(o, sv_type, r))
            plt.close()
            # Plotting total count over all samples
            fig = plt.figure(figsize =(10, 10))
            plt.bar(x_postion, sum, width=width, edgecolor='black')
            #plt.errorbar(x_postion, mean, yerr=stddev, fmt='none', ecolor='Black', elinewidth=2)
            plt.xlabel('Size (in bp)', fontsize=15)
            plt.ylabel('Total Count over all Samples', fontsize=15)
            plt.xticks(x_postion+0.5, x_labels, rotation=60, fontsize=5)
            plt.yticks(fontsize=10)
            plt.title(sv_type, fontsize=20)
            plt.tight_layout()
            plt.savefig('%s/total_count-size-distribution_%s_%s.png'%(o, sv_type, r))
            plt.close()
        # Plotting Stacked Bar Plots
        fig = plt.figure(figsize =(10, 10))
        plt.bar(x_postion, mean_list[0], width=width, edgecolor='black', color='grey', label='Insertions')
        plt.bar(x_postion, mean_list[1], width=width, bottom=mean_list[0], edgecolor='black', color='salmon', label='Deletions')
        plt.bar(x_postion, mean_list[2], width=width, bottom=list(np.array(mean_list[0]) + np.array(mean_list[1])), edgecolor='black', color='greenyellow', label='Others')
        #plt.errorbar(x_postion, mean, yerr=stddev, fmt='none', ecolor='Black', elinewidth=2)
        plt.xlabel('Size (in bp)', fontsize=15)
        plt.ylabel('Average Count per Sample', fontsize=15)
        plt.xticks(x_postion+0.5, x_labels, rotation=60, fontsize=5)
        plt.yticks(fontsize=10)
        plt.title("Stacked Bar Plot", fontsize=20)
        plt.legend()
        plt.tight_layout()
        plt.savefig('%s/stacked-average_count-size-distribution_%s.png'%(o, r))
        plt.close()
        fig = plt.figure(figsize =(10, 10))
        plt.bar(x_postion, sum_list[0], width=width, edgecolor='black', color='grey', label='Insertions')
        plt.bar(x_postion, sum_list[1], width=width, bottom=sum_list[0], edgecolor='black', color='salmon', label='Deletions')
        plt.bar(x_postion, sum_list[2], width=width, bottom=list(np.array(sum_list[0]) + np.array(sum_list[1])), edgecolor='black', color='greenyellow', label='Others')
        #plt.errorbar(x_postion, mean, yerr=stddev, fmt='none', ecolor='Black', elinewidth=2)
        plt.xlabel('Size (in bp)', fontsize=15)
        plt.ylabel('Total Count over all Samples', fontsize=15)
        plt.xticks(x_postion+0.5, x_labels, rotation=60, fontsize=5)
        plt.yticks(fontsize=10)
        plt.title("Stacked Bar Plot", fontsize=20)
        plt.legend()
        plt.tight_layout()
        plt.savefig('%s/stacked-total_count-size-distribution_%s.png'%(o, r))
        plt.close()
    
    binsize_1 = 10
    binsize_2 = 100
    binsize_3 = 1000
    binsize_4 = 10000
    n_bins_1 = 1000//binsize_1
    n_bins_2 = (10000-1000)//binsize_2
    n_bins_3 = (100000-10000)//binsize_3
    n_bins_4 = (1000000-100000)//binsize_4
    size_dist_1 = {'COMPLEX': [[] for _ in range(n_bins_1)], 'DEL': [[] for _ in range(n_bins_1)], 'INS': [[] for _ in range(n_bins_1)]}    # 0bp to 1,000bp
    size_dist_2 = {'COMPLEX': [[] for _ in range(n_bins_2)], 'DEL': [[] for _ in range(n_bins_2)], 'INS': [[] for _ in range(n_bins_2)]}    # 1,000bp to 10,000bp
    size_dist_3 = {'COMPLEX': [[] for _ in range(n_bins_3)], 'DEL': [[] for _ in range(n_bins_3)], 'INS': [[] for _ in range(n_bins_3)]}    # 10,000bp to 100,000bp
    size_dist_4 = {'COMPLEX': [[] for _ in range(n_bins_4)], 'DEL': [[] for _ in range(n_bins_4)], 'INS': [[] for _ in range(n_bins_4)]}    # 100,000bp to 1,000,000bp
    
    for sample in size_dist.keys():
        count_0 = 0
        count_1 = {'COMPLEX': [0 for _ in range(n_bins_1)], 'DEL': [0 for _ in range(n_bins_1)], 'INS': [0 for _ in range(n_bins_1)]}    # 0bp to 1,000bp
        count_2 = {'COMPLEX': [0 for _ in range(n_bins_2)], 'DEL': [0 for _ in range(n_bins_2)], 'INS': [0 for _ in range(n_bins_2)]}    # 1,000bp to 10,000bp
        count_3 = {'COMPLEX': [0 for _ in range(n_bins_3)], 'DEL': [0 for _ in range(n_bins_3)], 'INS': [0 for _ in range(n_bins_3)]}    # 10,000bp to 100,000bp
        count_4 = {'COMPLEX': [0 for _ in range(n_bins_4)], 'DEL': [0 for _ in range(n_bins_4)], 'INS': [0 for _ in range(n_bins_4)]}    # 100,000bp to 1,000,000bp
        data = size_dist[sample]
        for sv in data:
            t = sv.type
            if t == 'SNV':
                continue
            l = sv.len
            c = sv.count
            if l < 1000:
                bin = l//binsize_1
                count_1[t][bin] += c
                pass
            elif l < 10000:
                bin = (l-1000)//binsize_2
                count_2[t][bin] += c
                pass
            elif l < 100000:
                bin = (l-10000)//binsize_3
                count_3[t][bin] += c
                pass
            elif l < 1000000:
                bin = (l-100000)//binsize_4
                count_4[t][bin] += c
                pass
            else:
                count_0 += c
        #print("Found %d variants in sample %s of size more than 1Mbp"%(count_0, sample))
        for sv_type in ['COMPLEX', 'DEL', 'INS']:
            for i in range(n_bins_1):
                size_dist_1[sv_type][i].append(count_1[sv_type][i])
            for i in range(n_bins_2):
                size_dist_2[sv_type][i].append(count_2[sv_type][i])
            for i in range(n_bins_3):
                size_dist_3[sv_type][i].append(count_3[sv_type][i])
            for i in range(n_bins_4):
                size_dist_4[sv_type][i].append(count_4[sv_type][i])
    
    construct_plot(size_dist_1, '0-1kbp', output, binsize_1, 0)
    construct_plot(size_dist_2, '1kbp-10kbp', output, binsize_2, 1000)
    construct_plot(size_dist_3, '10kbp-100kbp', output, binsize_3, 10000)
    construct_plot(size_dist_4, '100kbp-1Mbp', output, binsize_4, 100000)


def plot_rausch_curve(allele_counts, plot_title, output, fname):
    filtered_allele_counts = [a.ac for a in allele_counts if a.ac != 0]
    #print("%d Variant Records filtered out due to AC = 0"%(len(allele_counts) - len(filtered_allele_counts)))
    counter = Counter(filtered_allele_counts)
    sorted_counter = sorted(counter.items(), key=lambda pair: pair[0])
    x = []
    y = []
    for i,j in sorted_counter:
        x.append(i)
        y.append(j)
    fig, ax = plt.subplots(figsize = (20,10))
    plt.plot(x, y)
    plt.xlabel("Variant Allele Count", fontsize=15)
    plt.ylabel("Number of Variant Sites", fontsize=15)
    plt.title(plot_title, fontsize=20)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig('%s/rausch-curve_%s.png'%(output, fname))
    plt.close()

SV = namedtuple('SV','name type bub len')
AC = namedtuple('AC','name type ac len')
SIZE_DIST = namedtuple('SIZED_DIST', 'type count len')

parser = argparse.ArgumentParser(prog='plot-qc-curves.py', description="Plotting QC curve from multisample biallelic vcf")
parser.add_argument('-vcf', metavar='VCF', help='multisample biallelic vcf')
parser.add_argument('-meta', metavar='META', help='sample metadata')
parser.add_argument('-output', metavar='OUTPUT', help='output directory')
args = parser.parse_args()

# initializing vcf
vcf = VCF(args.vcf)
vcf_samples = vcf.samples

# reading metadata (sample-to-population map and color coding)
metadata = pandas.read_csv(args.meta, sep='\t', header=0)
metadata = metadata[["Sample name", "Population code", "Superpopulation code"]]
metadata = metadata[metadata["Sample name"].isin(vcf_samples)].sort_values(by=["Superpopulation code"], ascending=False)

# sorted sample list required for Audano curve
non_AFR_samples = list(metadata[metadata["Superpopulation code"] != 'AFR']['Sample name'].values)
AFR_samples = list(metadata[metadata["Superpopulation code"] == 'AFR']['Sample name'].values)
random.shuffle(non_AFR_samples)
sorted_vcf_samples = non_AFR_samples + AFR_samples

# things to store
svs_per_sample = {}
ac = []
size_dist = {}

for sample in sorted_vcf_samples:
    svs_per_sample[sample] = []
    size_dist[sample] = []
count_unknown = 0
for variant in vcf:
    if variant.num_het+variant.num_hom_ref+variant.num_hom_alt == 0:
        count_unknown += 1
        continue
    name = variant.ID
    t, b, l = process_SV(name)
    sv = SV(name, t, b, int(l))
    allele_count = 0
    for index, g in enumerate(variant.genotypes):
        if g[0] != -1:
            allele_count += g[0] + g[1]
        s = vcf_samples[index]
        if g[0]+g[1] >= 1:
            svs_per_sample[s].append(sv)
            size_dist[s].append(SIZE_DIST(t, g[0]+g[1], int(l)))

    ac.append(AC(name, t, allele_count, int(l)))

# plotting the audano curve for the SVs
plot_audano_curve(svs_per_sample, metadata, args.output)

# plotting the rausch curve for different variant sets
plot_rausch_curve(ac, "All Variants with AC > 0", args.output, 'All-All')
plot_rausch_curve([a for a in ac if a.type == "INS"], "All INS with AC > 0", args.output, 'All-INS')
plot_rausch_curve([a for a in ac if a.type == "DEL"], "All DEL with AC > 0", args.output, 'All-DEL')
plot_rausch_curve([a for a in ac if a.type == "COMPLEX"], "All COMPLEX with AC > 0", args.output, 'All-COMPLEX')
plot_rausch_curve([a for a in ac if int(a.len) > 49], "All SVs with AC > 0", args.output, 'SV-All')
plot_rausch_curve([a for a in ac if (a.type == "INS" and int(a.len) > 49)], "All SV INS with AC > 0", args.output, 'SV-INS')
plot_rausch_curve([a for a in ac if (a.type == "DEL" and int(a.len) > 49)], "All SV DEL with AC > 0", args.output, 'SV-DEL')
plot_rausch_curve([a for a in ac if (a.type == "COMPLEX" and int(a.len) > 49)], "All SV COMPLEX with AC > 0", args.output, 'SV-COMPLEX')

# plotting the size distributions on different ranges
plot_sv_size_distribution(size_dist, args.output)