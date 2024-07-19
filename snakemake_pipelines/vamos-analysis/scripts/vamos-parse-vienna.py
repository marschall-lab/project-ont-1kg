import argparse
import sys
from collections import Counter

def run(vcf = None, sample = None, output = None):
    
    if output is None:
        writer = sys.stdout
    else:
        writer = open(output, 'w')

    sample_to_idx = []
    # vcf should not be gzipped
    with open(vcf, 'r') as vcffile:
        for line in vcffile:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                samples = line.rstrip().split('\t')[9:]
                for idx, s in enumerate(samples):
                    if sample != s:
                        continue
                    sample_to_idx.append(9+idx)
                continue
            line = line.rstrip().split('\t')
            info_field = line[7]
            rus = None
            altanno = []
            for info in info_field.split(';'):
                if info.startswith('RU='):
                    rus = info.split('=')[1].split(',')
                if info.startswith('ALTANNO_H1='):
                    assert rus != None
                    altanno.append(info.split('=')[1])
                if info.startswith('ALTANNO_H2='):
                    assert rus != None
                    altanno.append(info.split('=')[1])
            counter = Counter()
            for index in sample_to_idx:
                gts = line[index].split('/')
                for gt in gts:
                    counter.update([int(x) for x in altanno[int(gt)-1].split(',')])
            chrom = line[0]
            pos = line[1]
            for ru_idx, ru in enumerate(rus):
                if ru_idx in counter:
                    print(f"{chrom}\t{pos}\t{ru}\t{counter[ru_idx]}", file = writer)
                else:
                    print(f"{chrom}\t{pos}\t{ru}\t0", file = writer)
    
    if output != None:
        writer.close()

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='vamos-parse-vienna.py', description="Collects the repeating units from the ONT Vienna VCF and converts it into a TSV with counts")
    parser.add_argument("-vcf", required=True, help="Vamos VCF file")
    parser.add_argument("-sample", required=True, help="Name of the sample.")
    parser.add_argument("-output", help="Output TSV file.")

    options = parser.parse_args()

    run(**vars(options))