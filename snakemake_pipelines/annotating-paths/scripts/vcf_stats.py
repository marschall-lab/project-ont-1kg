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


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, help="VCF file")
    
    options = parser.parse_args()
    
    vcf = options.vcf

    if is_file_gzipped(vcf):
        reader = gzip.open(vcf, 'rt')
    else:
        reader = open(vcf, 'r')
    
    while True:
        line = reader.readline()
        if not line:
            break
        if line[0] == "#":
            continue
    
        line = line.split('\t')
        print(line[9:])
        for gt in line[9:]:
            pass
    

def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'

if __name__ == "__main__":
    main()