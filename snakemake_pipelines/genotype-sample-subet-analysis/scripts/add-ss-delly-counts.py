import argparse
import sys
import pandas as pd

def run(tsv = None, delly = None):
    
    tsv_reader = open(tsv, 'r')
    lines = []
    for line in tsv_reader:
        line = line.strip().split('\t')
        assert len(line) == 2
        lines.append(line)
    tsv_reader.close()
    delly_reader = open(delly, 'r')
    for line in delly_reader:
        line = line.strip().split('\t')
        assert len(line) == 1
        lines.append(['ss_delly', line])
    delly_reader.close()

    for line in lines:
        print(f'{line[0]}\t{line[1]}')

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='add-ss-delly-counts.py', description="adding sv counts from single sample delly calls")
    parser.add_argument("-tsv", required=True, help="Original output tsv")
    parser.add_argument("-delly", required=True, help="Delly output")

    options = parser.parse_args()

    run(**vars(options))