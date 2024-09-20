import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def run(bed=None, panel=None, callset=None, outdir=None):
    
    stats_writer = open(outdir+'stats.txt', 'w')
    
    # [<bub in gfa>, <bub has ps_hap allele in panel>, <bub in panel>]
    is_bub_present = defaultdict(lambda: [False, False, False])
    # [<sv is ps_hap allele in panel>, <sv in panel>, <sv in callset>]
    is_sv_present = defaultdict(lambda: [False, False, False])
    # reading the bed file and storing list of bubbles
    bed_reader = open(bed, 'r')
    count = 0
    for line in bed_reader:
        if line[0] == '#':
            continue
        _, _, _, bub_id = line.strip().split('\t')
        is_bub_present[bub_id][0] = True
        count += 1
    bed_reader.close()
    print(f'Number of bubbles in the rGFA = {count}', file=stats_writer)
    
    # reading the panel vcf
    samples=[]
    panel_reader = open(panel, 'r')
    for line in panel_reader:
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            samples=line.strip().split('\t')[9:]
            ps_start=0
            for index, sample in enumerate(samples):
                if 'path' in sample:
                    ps_start=9+index
                    break
            continue
        ps_gts=line.strip().split()[ps_start:]
        info_fields=line.split('\t')[7].split(';')
        for info in info_fields:
            if info.startswith('ID'):
                sv_id=info[3:]
                break
        bub_id=line.split()[2]
        if '1' in ps_gts:
            #assert is_bub_present[bub_id][0]==True
            is_bub_present[bub_id][1]=True
            is_sv_present[sv_id][0]=True
        is_bub_present[bub_id][2]=True
        is_sv_present[sv_id][1]=True
    panel_reader.close()

    # reading the callset vcf
    callset_reader = open(callset, 'r')
    for line in callset_reader:
        if line.startswith('#'):
            continue
        sv_id=line.split('\t')[2]
        if is_sv_present[sv_id][0]==False:
            continue
        is_sv_present[sv_id][2]=True
    callset_reader.close()

    count_panel = 0
    count_callset = 0
    sv_in_panel = {sv_type: [] for sv_type in ['SNV', 'INS', 'DEL', 'COMPLEX']}
    sv_in_callset = {sv_type: [] for sv_type in ['SNV', 'INS', 'DEL', 'COMPLEX']}
    for sv_id, x in is_sv_present.items():
        is_sv_in_panel, _, is_sv_in_callset = x
        _, _, sv_type, _, sv_len = sv_id.split('-')
        if is_sv_in_panel == False:
            continue
        count_panel += 1
        sv_in_panel[sv_type].append(int(sv_len))
        if is_sv_in_callset == True:
            count_callset += 1
            sv_in_callset[sv_type].append(int(sv_len))
    
    print(f'Number of pseudohaplotype SVs in the panel: {count_panel}', file=stats_writer)
    print(f'Number of the pseudohaplotype SVs that have been genotyped: {count_callset}', file=stats_writer)

    stats_writer.close()

    # Plotting counts of different sv types
    sv_types = ['INS', 'DEL', 'COMPLEX']
    sv_numbers = {'in_panel': [len([y for y in sv_in_panel[x] if y > 49]) for x in sv_types], 'in_callset': [len([y for y in sv_in_callset[x] if y > 49]) for x in sv_types]}
    x_pos = np.arange(len(sv_types))
    width = 0.4
    multiplier = 0

    fig, ax = plt.subplots(layout='constrained', figsize=(10,10))

    for x, y in sv_numbers.items():
        offset = width*multiplier
        rects = ax.bar(x_pos+offset, y, width, label=x)
        ax.bar_label(rects, padding=3)
        multiplier += 1
    ax.set_ylabel('Number of SVs')
    ax.set_title('Pseudohaplotype SV Count in Panel vs Callset')
    ax.set_xticks(x_pos+(width/2), sv_types)
    ax.legend(loc='upper right')
    ymin, ymax = ax.get_ylim()
    ax.set_ylim([0, ymax])
    plt.savefig(outdir+'counts.svg')
    
    # Plotting length histograms
    fig, ax = plt.subplots(1, 3, figsize=(30,10), tight_layout=True, sharey=True)
    
    for index, type in enumerate(['INS', 'DEL', 'COMPLEX']):
        ax[index].hist(sv_in_panel[type], bins=range(0, max(sv_in_panel[type])+1000, 1000), color='blue')
        ax[index].hist(sv_in_callset[type], bins=range(0, max(sv_in_callset[type])+1000, 1000), color='orange')
        ax[index].set_title(type)
    plt.yscale('log')
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='blue', lw=4), Line2D([0], [0], color='orange', lw=4)]
    plt.suptitle('SV Length Distribution for Pseudohaplotype SVs in Panel vs Callset')
    plt.legend(custom_lines, ['in_panel', 'in_callset'])
    plt.savefig(outdir+'length-dist.svg')

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='extract-bubbles.py', description="Extract bubbles from a tagged rGFA and outputs a BED file")
    parser.add_argument("-bed", required=True, help="Bubbles BED file")
    parser.add_argument("-panel", required=True, help="Biallelic Panel VCF")
    parser.add_argument("-callset", required=True, help="Biallelic Callset VCF")
    parser.add_argument("-outdir", required=True, help="Output directory")
    
    options = parser.parse_args()

    run(**vars(options))