'''
Produce stats for VCF created by pipeline.
'''

import sys
import argparse
import logging
from collections import defaultdict, namedtuple
import functools
import gzip
import re
import matplotlib.pyplot as plt


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, help="VCF file")
    parser.add_argument("--outdir", required=True, help="Output directory (direcotry should already exist)")
    
    options = parser.parse_args()
    
    vcf = options.vcf

    if is_file_gzipped(vcf):
        reader = gzip.open(vcf, 'rt')
    else:
        reader = open(vcf, 'r')
    
    gt_count = {}
    for c in range(45):
        gt_count[c] = 0
    sample_list = None
    sample_noGT_count = {}
    while True:
        line = reader.readline()
        if not line:
            break
        if line[0] == "##":
            continue
        if line[0] == "#":
            sample_list = line.rstrip().split('\t')[9:]
            for sample in sample_list:
                sample_noGT_count[sample] = 0
            continue

        line = line.rstrip().split('\t')
        count = 0
        for ind,gt in enumerate(line[9:]):
            gt = gt.split("|")
            if gt[0] == "." or gt[1] == ".":
                sample_noGT_count[sample_list[ind]] += 1
                continue
            else:
                count += 1
        gt_count[count] += 1
    
    plot_hist_numberOfAvailableSamples(gt_count, options.outdir)
    plot_dot_unknownBubblesPerSample(sample_noGT_count, options.outdir)


def plot_hist_numberOfAvailableSamples(gt_count, out):
    x = []
    y = []
    for n, count in gt_count.items():
        x.append(n)
        y.append(count)
    plt.bar(x,y)
    plt.yscale('log')
    plt.ylabel("Count of Bubbles")
    plt.xlabel("Number of Samples with Diploid data available for the Bubble")
    plt.savefig(out+"/hist_numberOfAvailableSamples.png")
    plt.close()
    

def plot_dot_unknownBubblesPerSample(sample_noGT_count, out):
    color_code = ["#ffd845", "#710027", "#778500", "#c44cfd"]
    
    # AFR = 0
    # AMR = 1
    # EAS = 2
    # SAS = 3

    pop = {'HG01123': 1,
           'HG03453': 0, 
           'HG02630': 0, 
           'HG02572': 0, 
           'HG01928': 1, 
           'HG03579': 0, 
           'HG00741': 1, 
           'HG02257': 0, 
           'HG00735': 1, 
           'HG00733': 1, 
           'HG00621': 2, 
           'HG02145': 0, 
           'HG03516': 0, 
           'HG01175': 1, 
           'HG01978': 1, 
           'HG03486': 0, 
           'HG02622': 0, 
           'HG03540': 0, 
           'HG01258': 1, 
           'HG03492': 3, 
           'HG02717': 0, 
           'HG02080': 2, 
           'HG02148': 1, 
           'HG01106': 1, 
           'HG01109': 1, 
           'HG01243': 1, 
           'NA18906': 0, 
           'HG01071': 1, 
           'HG00673': 2, 
           'HG02109': 0, 
           'HG01952': 1, 
           'HG02559': 0, 
           'HG01361': 1, 
           'HG03098': 0, 
           'HG02055': 0, 
           'HG02486': 0, 
           'HG02723': 0, 
           'HG02818': 0, 
           'NA21309': 0, 
           'HG01358': 1, 
           'HG02886': 0, 
           'HG01891': 0, 
           'HG00438': 2, 
           'NA20129': 0}
    
    color = []
    values = []
    pos = []
    p = 1
    for pop_code in [0,1,2,3]:
        for sample, value in sample_noGT_count.items():
            if pop[sample] != pop_code:
                continue
            color.append(color_code[pop_code])
            values.append(value)
            pos.append(p)
            p+=1

    plt.bar(pos, values, color=color)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.ylabel("Number of Unknown Bubbles")
    plt.xlabel("Samples")
    plt.savefig(out+"/hist_unknownBubblesPerSample.png")
    plt.close()


def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'

if __name__ == "__main__":
    main()