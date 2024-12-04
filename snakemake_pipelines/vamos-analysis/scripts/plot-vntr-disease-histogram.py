import argparse
import matplotlib.pyplot as plt

def plot_histogram(h, o, output, title):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    plt.hist(o, 20, density=True, facecolor='blue', alpha=0.5, label='ONT')
    plt.hist(h, 20, density=True, facecolor='green', alpha=0.5, label='HGSVC')
    plt.legend()
    plt.title(title)
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
    abca7_h = hgsvc_bps['ABCA7']
    abca7_o = ont_bps['ABCA7']
    plot_histogram(h=abca7_h, o=abca7_o, output=output+'abca7-histogram.svg', title='ABCA7 (in bps)')

    # plotting the PLIN4 RUs
    plin4_h = hgsvc_num_rus['PLIN4']
    plin4_o = ont_num_rus['PLIN4']
    plot_histogram(h=plin4_h, o=plin4_o, output=output+'plin4-histogram.svg', title='PLIN4 (in RUs)')


if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='plot-vntr-disease-histogram.py', description="Creates histogram for disease locus")
    parser.add_argument("-hgsvc", required=True, help="Stat files for HGSVC VNTRs")
    parser.add_argument("-ont", required=True, help="Stat files for ONT VNTRs")
    parser.add_argument("-output", required=True, help="Output directory")

    options = parser.parse_args()

    run(**vars(options))