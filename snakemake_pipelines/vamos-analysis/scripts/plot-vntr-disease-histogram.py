import argparse
import matplotlib.pyplot as plt

def plot_histogram(h, o, output, title, bins):
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(7, 7))
    # Plotting the RUs
    ax[0].hist(o[0], bins[0], density=True, facecolor='blue', alpha=0.5, label='ONT')
    ax[0].hist(h[0], bins[0], density=True, facecolor='green', alpha=0.5, label='HGSVC')
    ax[0].legend()
    ax[0].set_title(title+' (in RUs)')
    
    # Plotting the base pairs
    ax[1].hist(o[1], bins[1], density=True, facecolor='blue', alpha=0.5, label='ONT')
    ax[1].hist(h[1], bins[1], density=True, facecolor='green', alpha=0.5, label='HGSVC')
    ax[1].legend()
    ax[1].set_title(title+' (in base pairs)')

    fig.tight_layout()
    plt.savefig(output)
    pass

def check_vntr(sites, line):
    for name, pos in sites.items():
        if pos[0] == line[0] and pos[1] == line[1]:
            return name
    return None

def get_statistics(stats):
    disease_sites = {'ABCA7': ['chr19', '1012105'], 'PLIN4': ['chr19', '4494323']}
    num_rus = {x: [] for x in disease_sites.keys()}
    bps = {x: [] for x in disease_sites.keys()}
    for stat_file in stats.split(','):
        with open(stat_file, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                line = line.rstrip().split('\t')
                name=check_vntr(disease_sites, line)
                if name == None:
                    continue
                num_ru = [int(i) for i in line[2].split(',') if i != '.']
                bp = [int(i) for i in line[3].split(',') if i != '.']
                num_rus[name].extend(num_ru)
                bps[name].extend(bp)
    return num_rus, bps

def run(hgsvc, ont, output):
    hgsvc_num_rus, hgsvc_bps = get_statistics(hgsvc)
    ont_num_rus, ont_bps = get_statistics(ont)

    # plotting the ABCA7 basepairs
    bins_bps = [500*i for i in range(21)]
    bins_num_rus = [i for i in range(0, 401, 20)]
    plot_histogram(h=(hgsvc_num_rus['ABCA7'], hgsvc_bps['ABCA7']), o=(ont_num_rus['ABCA7'], ont_bps['ABCA7']), output=output+'abca7-histogram.svg', title='ABCA7', bins=(bins_num_rus, bins_bps))

    # plotting the PLIN4 RUs
    bins_num_rus = [i for i in range(10, 45, 2)]
    bins_bps = [i for i in range(1000, 4401, 200)]
    plot_histogram(h=(hgsvc_num_rus['PLIN4'], hgsvc_bps['PLIN4']), o=(ont_num_rus['PLIN4'], ont_bps['PLIN4']), output=output+'plin4-histogram.svg', title='PLIN4', bins=(bins_num_rus, bins_bps))


if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='plot-vntr-disease-histogram.py', description="Creates histogram for disease locus")
    parser.add_argument("-hgsvc", required=True, help="Stat files for HGSVC VNTRs")
    parser.add_argument("-ont", required=True, help="Stat files for ONT VNTRs")
    parser.add_argument("-output", required=True, help="Output directory")

    options = parser.parse_args()

    run(**vars(options))