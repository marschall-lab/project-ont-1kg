import argparse
import pandas
import sys
import matplotlib.pyplot as plt
import matplotlib
from pysam import VariantFile

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='sv_count.py', description="Count SVs from list of VCFs")
    parser.add_argument('-meta', metavar='META', help='metadata')
    parser.add_argument('-vcf', metavar='VCF', help='comma-separated list of vcfs')
    parser.add_argument('-output', metavar='OUTPUT', help='output directory')
    args = parser.parse_args()
    
    # reading metadata (sample-to-population map and color coding)
    metadata = pandas.read_csv(args.meta, sep='\t', header=0)
    metadata = metadata[["Sample name", "Population code", "Superpopulation code"]]
    
    # reading vcfs and extracting genotypes
    vcfs = args.vcf.split(',')
    data = {}
    for n,v in enumerate(vcfs):
        sample=v[-14:-7]
        sample_data=metadata[metadata["Sample name"] == sample]
        d = [0, 0, sample_data["Population code"].values[0], sample_data["Superpopulation code"].values[0]]      # first element counts number of variant positions with SVs. second element counts SVs.
        reader = VariantFile(v)
        assert len(reader.header.samples) == 1
        assert reader.header.samples[0] == sample
        for rec in reader.fetch():
            gts = [s['GT'] for s in rec.samples.values()]
            assert len(gts) == 1
            gt = gts[0]
            if gt[0] == None:
                continue
            if gt[0] == 0  and gt[1] == 0:
                continue
            elif gt[0] != 0 and gt[1] != 0:
                d[0] += 1
                d[1] += 2
            else:
                d[0] += 1
                d[1] += 1
        data[sample] = d
        reader.close()
        if (n+1)%10==0:
            print("Read %d VCFs..."%(n+1), file=sys.stderr)
    
    # sort data by them superpopulation code
    sorted_data = dict(sorted(data.items(), key=lambda x:x[1][3]))
    sv_count_by_pop = {'AFR': [], 'AMR': [], 'EAS': [], 'EUR': [], 'SAS': []}
    var_count_by_pop = {'AFR': [], 'AMR': [], 'EAS': [], 'EUR': [], 'SAS': []}
    color_by_pop = {'AFR': '#ffd845', 'AMR': '#710027', 'EAS': '#778500', 'EUR': '#018ead', 'SAS': '#c44cfd'}
    sample_list = []
    color_list = []
    sv_count_list = []
    var_count_list = []
    for sample, value in sorted_data.items():
        sample_list.append(sample)
        color_list.append(color_by_pop[value[3]])
        sv_count_list.append(value[1])
        var_count_list.append(value[0])
        sv_count_by_pop[value[3]].append(value[1])
        var_count_by_pop[value[3]].append(value[0])
    
    # create legend
    handles = []
    labels = []
    color_index = 0
    for pop,color in color_by_pop.items():
        line = matplotlib.patches.Patch(color=color, label=pop)
        handles.append(line)
        labels.append(pop)
        color_index += 1
    
    # Plotting SV count sample-wise
    fig, ax = plt.subplots(figsize=(20,5))
    ax.bar(sample_list, sv_count_list, color=color_list)
    ax.set_ylabel('SV Count')
    ax.set_title('SV Counts Sample-Wise')
    plt.figlegend(handles, labels, framealpha=1, frameon=True)
    plt.tight_layout()
    plt.savefig(args.output+'/sv_count_sample-wise.png')
    
    # Plotting Variant Count sample-wise
    fig, ax = plt.subplots(figsize=(20,5))
    ax.bar(sample_list, var_count_list, color=color_list)
    ax.set_ylabel('Number of Variants')
    ax.set_title('Number of Variant Positions with an SV')
    plt.figlegend(handles, labels, framealpha=1, frameon=True)
    plt.tight_layout()
    plt.savefig(args.output+'/var_count_sample-wise.png')
    
    # Plotting SV count population-wise
    fig = plt.figure(figsize =(10, 10))
    ax = fig.add_subplot(111)
    bp = ax.boxplot([sv_count_by_pop[x] for x in sv_count_by_pop.keys()], patch_artist = True)
    colors = [v for _,v in color_by_pop.items()]
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    for whisker in bp['whiskers']:
        whisker.set(color ='black', linewidth = 1, linestyle =":")
    for cap in bp['caps']:
        cap.set(color ='black', linewidth = 1)
    for median in bp['medians']:
        median.set(color ='black', linewidth = 1)
    for flier in bp['fliers']:
        flier.set(marker ='D', color ='black', alpha = 0.5)    
    ax.set_xticklabels([x for x,_ in color_by_pop.items()])
    ax.set_ylabel('SV Count')
    ax.set_title('SV Counts Population-Wise')
    plt.figlegend(handles, labels, framealpha=1, frameon=True)
    plt.tight_layout()
    plt.savefig(args.output+'/sv_count_population-wise.png')
    
    # Plotting Variant count population-wise
    fig = plt.figure(figsize =(10, 10))
    ax = fig.add_subplot(111)
    bp = ax.boxplot([var_count_by_pop[x] for x in sv_count_by_pop.keys()], patch_artist = True)
    colors = [v for _,v in color_by_pop.items()]
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    for whisker in bp['whiskers']:
        whisker.set(color='black', linewidth=1, linestyle=":")
    for cap in bp['caps']:
        cap.set(color='black', linewidth=1)
    for median in bp['medians']:
        median.set(color='black', linewidth=1)
    for flier in bp['fliers']:
        flier.set(marker='D', color='black', alpha=0.5)    
    ax.set_xticklabels([x for x,_ in color_by_pop.items()])
    ax.set_ylabel('Number of Variants')
    ax.set_title('Number of Variant Positions with an SV')
    plt.figlegend(handles, labels, framealpha=1, frameon=True)
    plt.tight_layout()
    plt.savefig(args.output+'/var_count_population-wise.png')

    # Writing data into file
    writer = open(args.output+'/sv_count.tsv', 'w')
    print("Sample\tPopulation\tSuperpopulation\tSV Count\tNumber of Variant Positions with SVs", file=writer)
    for sample, value in sorted_data.items():
        print("%s\t%s\t%s\t%d\t%d"%(sample, value[2], value[3], value[1], value[0]), file=writer)
    writer.close()