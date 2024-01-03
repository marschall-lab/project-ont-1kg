#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(prog='genotyping-evaluation.py', description=__doc__)
parser.add_argument('directory', metavar='DIRECTORY', help='directory to the concordance summaries.')
parser.add_argument('regions', metavar='REGIONS', help='comma-separated list of regions (can be multiallelic or biallelic).')
parser.add_argument('vartype', metavar='VARTYPE', help='comma-separated list of variant types.')

args = parser.parse_args()
regions = args.regions.split(',')
regions.sort()
vartype = args.vartype.split(',')
vartype.sort()

title_row = ['variant_class','weighted_concordance', 'genotype_concordance_all', 'genotype_concordance_absent', 'genotype_concordance_het', 'genotype_concordance_hom', 'typed', 'total_baseline_variants', 'total_unmatched', 'precision', 'recall']
print('\t'.join(title_row))
for r in regions:
    for v in vartype:
        file = '%s/%s_%s/summary.txt'%(args.directory, r, v)
        result = ""
        result += "%s_%s"%(r,v)
        with open(file, 'r') as f:
            for n, line in enumerate(f):
                if n < 16:
                    continue
                val = float(line.rstrip().split('\t')[1])
                val = str(round(val, 4))
                result += "\t%s"%(val)
        print(result)
        
            

