import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def run(vcf=None, outdir=None):
    
    sv_lens = []
    all_sv_lens = []
    # reading the callset vcf
    callset_reader = open(vcf, 'r')
    for line in callset_reader:
        if line.startswith('#'):
            continue
        sv_id=line.split('\t')[2]
        aa_field=line.split('\t')[7].split(';')[-1]
        assert aa_field.startswith('ANCESTRAL_ALLELE')
        all_sv_lens.append(int(sv_id.split('-')[-1]))
        if aa_field.split('=')[-1] == '.':
            continue
        sv_lens.append(int(sv_id.split('-')[-1]))
    callset_reader.close()
    
    # Plotting line plot and histogram in same plot
    fig, ax = plt.subplots(1, 1, figsize=(10,10), tight_layout=True)
    n, x = np.histogram(all_sv_lens, bins=range(0, max(sv_lens)+1000, 1000))
    bin_centers = 0.5*(x[1:]+x[:-1])
    #plt.plot(bin_centers, n, color='blue', label='All SVs')
    ax.hist(all_sv_lens, bins=range(0, max(sv_lens)+1000, 1000), alpha=0.5 , color='orange', label='All SVs')
    n, x = np.histogram(sv_lens, bins=range(0, max(sv_lens)+1000, 1000))
    bin_centers = 0.5*(x[1:]+x[:-1])
    plt.plot(bin_centers, n, color='black', label='SVs with known Ancestral All') ## using bin_centers rather than edges
    ax.set_title('Lengths of SVs annotated with an Ancestral Allele', fontsize=16)
    plt.xlabel('Lengths (in bp)', fontsize=18)
    plt.ylabel('Counts', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(prop={'size': 16})
    plt.yscale('log')
    plt.savefig(outdir+'sv-length-dist.svg')

    """
    # Plotting length histograms
    fig, ax = plt.subplots(1, 1, figsize=(10,10), tight_layout=True)
    ax.hist(all_sv_lens, bins=range(0, max(sv_lens)+1000, 1000), color='blue', label='All SVs')
    ax.hist(sv_lens, bins=range(0, max(sv_lens)+1000, 1000), color='orange', label='SVs with known Ancestral All')
    ax.set_title('Lengths of SVs annotated with an Ancestral Allele', fontsize=16)
    plt.xlabel('Lengths', fontsize=18)
    plt.ylabel('Counts', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(prop={'size': 16})
    plt.savefig(outdir+'sv-length-dist.svg')
    
    # Plotting length line histograms
    fig, ax = plt.subplots(1, 1, figsize=(10,10), tight_layout=True)
    n, x = np.histogram(all_sv_lens, bins=range(0, max(sv_lens)+1000, 1000))
    bin_centers = 0.5*(x[1:]+x[:-1])
    #plt.plot(bin_centers, n, color='blue', label='All SVs')
    ax.hist(all_sv_lens, bins=range(0, max(sv_lens)+1000, 1000), color='blue', label='All SVs')
    n, x = np.histogram(sv_lens, bins=range(0, max(sv_lens)+1000, 1000))
    bin_centers = 0.5*(x[1:]+x[:-1])
    plt.plot(bin_centers, n, color='orange', label='SVs with known Ancestral All') ## using bin_centers rather than edges
    ax.set_title('Lengths of SVs annotated with an Ancestral Allele', fontsize=16)
    plt.xlabel('Lengths', fontsize=18)
    plt.ylabel('Counts', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(prop={'size': 16})
    plt.yscale('log')
    plt.savefig(outdir+'sv-length-dist-line.svg')
    """

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='plot-sv-lengths.py', description="Plot the SV lengths of the VCF records annotated by the pipeline")
    parser.add_argument("-vcf", required=True, help="Annotated VCF")
    parser.add_argument("-outdir", required=True, help="Output directory")
    
    options = parser.parse_args()

    run(**vars(options))