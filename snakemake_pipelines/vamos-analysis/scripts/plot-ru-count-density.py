import argparse
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.stats import pearsonr
import seaborn

def run(miller = None, vienna = None, output = None):
    
    miller_reader = open(miller, 'r')
    vienna_reader = open(vienna, 'r')
    miller_ru = defaultdict(lambda: 0)
    vienna_ru = defaultdict(lambda: 0)
    while True:
        line = miller_reader.readline()
        if not line:
            break
        chrom, pos, ru, count = line.rstrip().split('\t')
        miller_ru[(chrom, int(pos), ru)] = int(count)
    while True:
        line = vienna_reader.readline()
        if not line:
            break
        chrom, pos, ru, count = line.rstrip().split('\t')
        vienna_ru[(chrom, int(pos), ru)] = int(count)
    miller_counts = []
    vienna_counts = []
    for info in vienna_ru.keys():
        if miller_ru[info] == vienna_ru[info] and vienna_ru[info] == 0:
            continue
        if miller_ru[info] > 100 or vienna_ru[info] > 100:
            continue
        miller_counts.append(miller_ru[info])
        vienna_counts.append(vienna_ru[info])
    
    g = seaborn.JointGrid(x=miller_counts, y=vienna_counts, height=10, ratio=5, ylim=(0, 100), xlim=(0,100))
    g.plot_marginals(seaborn.histplot, bins=10)
    g.plot_joint(plt.hexbin, bins='log', gridsize=50, cmap='Blues')
    g.set_axis_labels("RU Counts from Miller VCF", "RU Counts from Vienna VCF")
    cbar_ax = g.fig.add_axes([.85, .2, .02, .6])  # x, y, width, height
    plt.colorbar(cax=cbar_ax)
    plt.title(f'PCC: {pearsonr(miller_counts, vienna_counts)}')
    g.savefig(output)
    
    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='plot-ru-count-density.py', description="Creates a scatter plot of the repeating unit counts comparing the vamos run on miller vcf and vienna vcf")
    parser.add_argument("-miller", required=True, help="Miller Vamos VCF file")
    parser.add_argument("-vienna", required=True, help="Vienna Vamos VCF file")
    parser.add_argument("-output", required=True, help="Output scatter plot.")

    options = parser.parse_args()

    run(**vars(options))