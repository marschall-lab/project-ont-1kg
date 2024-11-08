import argparse
import matplotlib.pyplot as plt
import matplotlib
import seaborn

def run(same_sample, diff_sample, output):
    
    same_sample_pccs = []
    for file in same_sample.split(','):
        with open(file, 'r') as reader:
            for line in reader:
                # labelled as PCC_RU-Scatter_Whole
                if line.startswith('01'):
                    same_sample_pccs.append(float(line.split('\t')[2]))         
    diff_sample_pccs = []
    for file in diff_sample.split(','):
        with open(file, 'r') as reader:
            for line in reader:
                # labelled as PCC_RU-Scatter_Whole
                if line.startswith('01'):
                    diff_sample_pccs.append(float(line.split('\t')[2]))
    
    #colors = seaborn.color_palette("tab10", 2)
    labels = ['Same-sample comparison', 'Diff-sample comparison']
    #handles = [matplotlib.patches.Patch(color=colors[i], label=[labels[i]]) for i in [0,1]]
    data = [same_sample_pccs, diff_sample_pccs]
    plt.boxplot(data, positions=[1, 3], labels=labels, notch=True, widths=1)
    #plt.boxplot(same_sample_pccs, positions=[1], patch_artist = True, notch=True, widths=1)
    #plt.boxplot(diff_sample_pccs, positions=[3], patch_artist = True, notch=True, widths=1)
    #plt.legend(handles, labels)
    #fig, ax = plt.subplots()
    plt.savefig(output)

    
    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='plot-hgsvc-qc.py', description="Creates plots for HGSVC-based QC")
    parser.add_argument("-same-sample", required=True, help="Comma-separated list of same-sample comparison stat files")
    parser.add_argument("-diff-sample", required=True, help="Comma-separated list of diff-sample comparison stat files")
    parser.add_argument("-output", required=True, help="Output boxplot")

    options = parser.parse_args()

    run(**vars(options))