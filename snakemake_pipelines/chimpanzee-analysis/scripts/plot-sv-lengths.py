import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def run(vcf=None, outdir=None):
    
    sv_lens = []
    # reading the callset vcf
    callset_reader = open(vcf, 'r')
    for line in callset_reader:
        if line.startswith('#'):
            continue
        sv_id=line.split('\t')[2]
        aa_field=line.split('\t')[7].split(';')[-1]
        assert aa_field.startswith('ANCESTRAL_ALLELE')
        if aa_field.split('=')[-1] == '.':
            continue
        sv_lens.append(int(sv_id.split('-')[-1]))
    callset_reader.close()
    
    # Plotting length histograms
    fig, ax = plt.subplots(1, 1, figsize=(10,10), tight_layout=True)
    ax.hist(sv_lens, bins=range(0, max(sv_lens)+1000, 1000), color='blue')
    ax.set_title('Lengths of SVs annotated with an ancestral allele')
    plt.yscale('log')
    plt.savefig(outdir+'sv-length-dist.svg')

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='plot-sv-lengths.py', description="Plot the SV lengths of the VCF records annotated by the pipeline")
    parser.add_argument("-vcf", required=True, help="Annotated VCF")
    parser.add_argument("-outdir", required=True, help="Output directory")
    
    options = parser.parse_args()

    run(**vars(options))