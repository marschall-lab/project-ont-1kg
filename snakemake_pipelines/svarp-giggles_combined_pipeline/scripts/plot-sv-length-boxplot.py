'''
This script plots similar to sv_count.py but instead of count, it plots the total sv lenghts found in the sample.
'''

import argparse
import pandas
import sys
import matplotlib.pyplot as plt
import matplotlib
from cyvcf2 import VCF

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='sv_count.py', description="Count SVs from list of VCFs")
    parser.add_argument('-meta', metavar='META', help='metadata')
    parser.add_argument('-vcf', metavar='VCF', help='multisample biallelic vcf')
    parser.add_argument('-table', metavar='TABLE', help='tsv file of data')
    parser.add_argument('-output', metavar='OUTPUT', help='output directory')
    args = parser.parse_args()
    
    # reading metadata (sample-to-population map and color coding)
    metadata = pandas.read_csv(args.meta, sep='\t', header=0)
    metadata = metadata[["Sample name", "Population code", "Superpopulation code"]]
    
    # reading vcfs and extracting genotypes
    if args.vcf:
        # initializing vcf
        vcf = VCF(args.vcf)
        vcf_samples = vcf.samples
        data = {}
        for s in vcf_samples:
            sample_data=metadata[metadata["Sample name"] == s]
            data[s] = [[0,0,0], [0,0,0], sample_data["Population code"].values[0], sample_data["Superpopulation code"].values[0]]      # first element counts number of HET SVs [COMPLEX, DEL, INS]. second element counts HOM SVs [COMPLEX, DEL, INS].  
            
        sv_type_to_index = {'COMPLEX': 0, 'DEL': 1, 'INS': 2}
        # Writing data into file
        writer = open(args.output+'/plot-sv-length-boxplot.tsv', 'w')
        print("Sample\tPopulation\tSuperpopulation\tHET_COMPLEX\tHET_DEL\tHET_INS\tHOM_COMPLEX\tHOM_DEL\tHOM_INS", file=writer)
        for variant in vcf:
            if variant.num_het+variant.num_hom_ref+variant.num_hom_alt == 0:
                continue
            _, _, sv_type, bub, length = variant.ID.split('-')
            length = int(length)
            if length < 50:
                continue
            for index, g in enumerate(variant.genotypes):
                sample = vcf_samples[index]
                if g[0] == -1:
                    continue
                assert g[0] in [0, 1] and g[1] in [0, 1]
                if g[0] == 0 and g[1] == 0:
                    continue
                if g[0]+g[1] == 1:
                    # HET SV    
                    data[sample][0][sv_type_to_index[sv_type]] += length
                elif g[0]+g[1] == 2:
                    # HOM SV
                    data[sample][1][sv_type_to_index[sv_type]] += length
        vcf.close()
        for sample in vcf_samples:
            d = data[sample]
            print("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d"%(sample, d[2], d[3], *d[0], *d[1]), file=writer)
        writer.close()
    
    elif args.table:
        data = {}
        reader = open(args.table, 'r')
        for line in reader:
            if line[0] == "S":
                continue
            sample, pop, spop, het_com, het_del, het_ins, hom_com, hom_del, hom_ins = line.rstrip().split('\t')
            data[sample] = [[int(het_com), int(het_del), int(het_ins)], [int(hom_com), int(hom_del), int(hom_ins)], pop, spop]
    else:
        sys.exit("No table or list of vcfs given")

    # sort data by them superpopulation code
    sorted_data = dict(sorted(data.items(), key=lambda x:x[1][3]))
    hom_count_by_pop = {'AFR': [], 'AMR': [], 'EAS': [], 'EUR': [], 'SAS': []}
    het_count_by_pop = {'AFR': [], 'AMR': [], 'EAS': [], 'EUR': [], 'SAS': []}
    color_by_pop = {'AFR': '#ffd845', 'AMR': '#710027', 'EAS': '#778500', 'EUR': '#018ead', 'SAS': '#c44cfd'}
    sample_list = []
    color_list = []
    hom_count_list = []
    het_count_list = []
    for sample, value in sorted_data.items():
        sample_list.append(sample)
        color_list.append(color_by_pop[value[3]])
        hom_count_list.append(sum(value[1])/1000000)
        het_count_list.append(sum(value[0])/1000000)
        hom_count_by_pop[value[3]].append(sum(value[1])/1000000)
        het_count_by_pop[value[3]].append(sum(value[0])/1000000)
    
    # create legend
    handles = []
    labels = []
    color_index = 0
    for pop,color in color_by_pop.items():
        line = matplotlib.patches.Patch(color=color, label=pop)
        handles.append(line)
        labels.append(pop)
        color_index += 1
    
    y_range = (0, 13)

    # Plotting SV count sample-wise
    fig, ax = plt.subplots(figsize=(10,10), dpi=100)
    ax.bar(sample_list, hom_count_list, color=color_list)
    ax.set_ylabel('HOM SV Length (in Mbp)')
    ax.set_title('HOM SV Counts Sample-Wise')
    plt.ylim(y_range)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.figlegend(handles, labels, framealpha=1, frameon=True)
    plt.tight_layout()
    plt.savefig(args.output+'/hom_length_sample-wise.png')
    
    # Plotting Variant Count sample-wise
    fig, ax = plt.subplots(figsize=(10,10), dpi=100)
    ax.bar(sample_list, het_count_list, color=color_list)
    ax.set_ylabel('HET SV Length (in Mbp)')
    ax.set_title('HET SV Counts Sample-Wise')
    plt.ylim(y_range)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.figlegend(handles, labels, framealpha=1, frameon=True)
    plt.tight_layout()
    plt.savefig(args.output+'/het_length_sample-wise.png')
    
    # Plotting SV count population-wise
    fig = plt.figure(figsize =(10, 10))
    ax = fig.add_subplot(111)
    bp = ax.boxplot([hom_count_by_pop[x] for x in hom_count_by_pop.keys()], patch_artist = True, notch=True)
    colors = [v for _,v in color_by_pop.items()]
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    for whisker in bp['whiskers']:
        whisker.set(color ='black', linewidth = 1)
    for cap in bp['caps']:
        cap.set(color ='black', linewidth = 1)
    for median in bp['medians']:
        median.set(color ='black', linewidth = 1)
    for flier in bp['fliers']:
        flier.set(marker ='o', color ='black', alpha = 1)    
    ax.set_xticklabels([x for x,_ in color_by_pop.items()])
    ax.set_ylabel('HOM SV Length (in Mbp)', fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.ylim(y_range)
    plt.tight_layout()
    plt.savefig(args.output+'/hom_length_population-wise.png')
    
    # Plotting Variant count population-wise
    fig = plt.figure(figsize =(10, 10))
    ax = fig.add_subplot(111)
    bp = ax.boxplot([het_count_by_pop[x] for x in het_count_by_pop.keys()], patch_artist = True, notch=True)
    colors = [v for _,v in color_by_pop.items()]
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    for whisker in bp['whiskers']:
        whisker.set(color ='black', linewidth = 1)
    for cap in bp['caps']:
        cap.set(color ='black', linewidth = 1)
    for median in bp['medians']:
        median.set(color ='black', linewidth = 1)
    for flier in bp['fliers']:
        flier.set(marker ='o', color ='black', alpha = 1)    
    ax.set_xticklabels([x for x,_ in color_by_pop.items()])
    ax.set_ylabel('HET SV Length (in Mbp)', fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.ylim(y_range)
    plt.tight_layout()
    plt.savefig(args.output+'/het_length_population-wise.png')