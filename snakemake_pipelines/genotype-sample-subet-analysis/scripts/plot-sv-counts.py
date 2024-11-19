import argparse
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

def read_tsv(file):
    sv_counter = {}
    with open(file, 'r') as f:
        for line in f:
            line=line.strip().split('\t')
            assert len(line) == 2
            sv_counter[line[0]] = [int(i) for i in line[1].split(',')]
    return sv_counter

def run(tsvs=None, output=None, title=None):
    tsvs = tsvs.split(',')
    ranges = [i.split('/')[-1].split('.')[0] for i in tsvs]
    ranges.sort(key = lambda range: int(range.split('-')[0]))
    #categories = ['callset', 'sniffles', 'delly', 'svarp', 'single_sample_delly']
    categories = ['single_sample_delly_AFR', 'delly_AFR', 'callset_AFR', 'single_sample_delly_non_AFR', 'delly_non_AFR', 'callset_non_AFR']
    # create legend
    handles = []
    labels = []
    colors = sns.color_palette("tab10", len(categories))
    count = 0
    for name in categories:
        line = matplotlib.patches.Patch(color=colors[count], label=name)
        handles.append(line)
        labels.append(name)
        count += 1
    # preparing data
    data = {}
    for file, cov_range in zip(tsvs, ranges):
        counter = read_tsv(file)
        data[cov_range] = counter
    # plotting
    fig = plt.figure(figsize =(20, 10))
    bp = {}
    tick_pos = [(i*(len(categories)+1))+(len(categories)+1)/2 for i in range(len(ranges))]
    for index, cat in enumerate(categories):
        d = [data[i][cat] for i in ranges]
        bp[cat] = plt.boxplot(d, positions=[(index+1)+(i*(len(categories)+1)) for i in range(len(ranges))], patch_artist=True, notch=False, widths=0.4)
        # styling the boxes
        for patch in bp[cat]['boxes']:
            patch.set_facecolor(colors[index])
        for whisker in bp[cat]['whiskers']:
            whisker.set(color ='black', linewidth = 1)
        for cap in bp[cat]['caps']:
            cap.set(color ='black', linewidth = 1)
        for median in bp[cat]['medians']:
            median.set(color ='black', linewidth = 1)
        for flier in bp[cat]['fliers']:
            flier.set(marker ='o', color ='black', alpha = 1)
    for index in range(len(ranges)):
        if index == 0:
            continue
        plt.axvline(x = index*(len(categories)+1), color = 'black', linewidth = 1, linestyle='--')
    
    plt.title(f'SV Count for Discovery vs Final Genotypes - {title}', fontsize=20)
    plt.xticks(tick_pos, ranges, fontsize=15)
    xmin, xmax, _, _ = plt.axis()
    plt.xlim([xmin-1, xmax+1])
    plt.xlabel(f'{title}', fontsize=20)
    plt.ylabel('Number of SVs', fontsize=20)
    plt.yticks(fontsize=15)
    plt.tight_layout()
    plt.legend(handles, labels, fontsize=15)
    plt.savefig(output)
    plt.close()


if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='count-svs-per-sample.py', description="count svs per sample")
    parser.add_argument("-tsvs", required=True, help="comma-separated list of tsvs with the count information")
    parser.add_argument("-output", required=True, help="output plot")
    parser.add_argument("-title", required=True, help="plot title")

    options = parser.parse_args()

    run(**vars(options))