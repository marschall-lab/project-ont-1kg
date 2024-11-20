import argparse
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import seaborn

def read_summary_file(file, c1, c2):
    data = {}
    for line in file:
        if line[0] == '#':
            line = line.strip().split('\t')
            i1 = line.index(c1)
            i2 = line.index(c2)
            continue
        line = line.strip().split('\t')
        chrom = line[0]
        pos = line[1]
        key = chrom+'_'+pos
        data[key] = float(line[i2])-float(line[i1])
    
    return data

def run(hgsvc, ont, range, output):
    
    title_map = {
        'max': 'Min - Max',
        '95': '5 Percentile - 95 Percentile',
        '75': '25 Percentile - 75 Percentile'
    }

    column_map = {
        'max': ['MIN_RUS', 'MAX_RUS'],
        '95': ['5%_RUS', '95%_RUS'],
        '75': ['25%_RUS', '75%_RUS']
    }

    title = title_map[range]

    col1 = column_map[range][0]
    col2 = column_map[range][1]

    hgsvc_data = None
    ont_data = None
    with open(hgsvc, 'r') as file:
        hgsvc_data = read_summary_file(file, col1, col2)
    with open(ont, 'r') as file:
        ont_data = read_summary_file(file, col1, col2)

    hgsvc_sites = list(hgsvc_data.keys())
    ont_sites = list(ont_data.keys())

    print(f'Number of hgsvc sites: {len(hgsvc_sites)}')
    print(f'Number of ont sites: {len(ont_sites)}')

    intersect_keys = list(set(hgsvc_sites).intersection(ont_sites))
    print(f'Number of intersecting sites: {len(intersect_keys)}')

    hgsvc_values = []
    ont_values = []
    hgsvc_values_subset100 = []
    ont_values_subset100 = []
    for key in intersect_keys:
        h = hgsvc_data[key]
        o = ont_data[key]
        hgsvc_values.append(h)
        ont_values.append(o)
        if h < 100 and o < 100:
            hgsvc_values_subset100.append(h)
            ont_values_subset100.append(o)
        
    # plotting the density of RU spread by taking the difference of different ranges

    print(f'Number of positions plotted in full plot: {len(hgsvc_values)}')
    g = seaborn.JointGrid(x=hgsvc_values, y=ont_values, height=10, ratio=5)
    g.plot_marginals(seaborn.histplot, bins=100)
    g.plot_joint(plt.hexbin, bins='log', gridsize=50, cmap='Blues')
    g.set_axis_labels("Range Difference in HGSVC", "Range Difference in ONT", fontsize=20)
    g.fig.suptitle(f'{title} - Full\nPCC: {pearsonr(hgsvc_values, ont_values)[0]}', fontsize=20)
    cbar_ax = g.fig.add_axes([.95, .2, .02, .6])  # x, y, width, height
    plt.colorbar(cax=cbar_ax)
    g.savefig(f'{output}-full.svg')
    
    print(f'Number of positions plotted in subset100 plot: {len(hgsvc_values_subset100)}')
    f = seaborn.JointGrid(x=hgsvc_values_subset100, y=ont_values_subset100, height=10, ratio=5, ylim=(-100, 100), xlim=(-100,100))
    #g = seaborn.JointGrid(x=hgsvc_counts, y=ont_counts, height=10, ratio=5)
    f.plot_marginals(seaborn.histplot, bins=range(0, 100))
    f.plot_joint(plt.hexbin, bins='log', gridsize=50, cmap='Blues')
    f.set_axis_labels("Range Difference in HGSVC", "Range Difference in ONT", fontsize=20)
    f.fig.suptitle(f'{title} - Subset100\nPCC: {pearsonr(hgsvc_values_subset100, ont_values_subset100)[0]}', fontsize=20)
    cbar_ax = f.fig.add_axes([.95, .2, .02, .6])  # x, y, width, height
    plt.colorbar(cax=cbar_ax)
    f.savefig(f'{output}-subset100.svg')
    
    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='plot-hgsvc-qc.py', description="Creates plots for HGSVC-based QC")
    parser.add_argument("-hgsvc", required=True, help="Summary statistics of HGSVC VNTRs")
    parser.add_argument("-ont", required=True, help="Summary statistics of ONT VNTRs")
    parser.add_argument("-range", required=True, help="Determines the range of values to use (max | 95 | 75)")
    parser.add_argument("-output", required=True, help="Output prefix")

    options = parser.parse_args()

    run(**vars(options))