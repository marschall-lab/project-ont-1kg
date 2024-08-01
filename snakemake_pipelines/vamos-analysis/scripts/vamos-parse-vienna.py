import argparse
from collections import Counter

def run(vcf = None, sample = None, output = None):
    
    count_writer = open(output+'-count.tsv', 'w')
    seq_writer = open(output+'-sequence.tsv', 'w')

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
            chrom = line[0]
            pos = line[1]
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
                    seq=''.join([rus[int(x)] for x in altanno[int(gt)-1].split(',')])
                    print(f"{chrom}\t{pos}\t{seq}", file=seq_writer)
                    counter.update([int(x) for x in altanno[int(gt)-1].split(',')])
            for ru_idx, ru in enumerate(rus):
                if ru_idx in counter:
                    print(f"{chrom}\t{pos}\t{ru}\t{counter[ru_idx]}", file = count_writer)
                else:
                    print(f"{chrom}\t{pos}\t{ru}\t0", file = count_writer)
    
    count_writer.close()
    seq_writer.close()

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='vamos-parse-vienna.py', description="Collects the repeating units from the ONT Vienna VCF and converts it into a TSV with counts")
    parser.add_argument("-vcf", required=True, help="Vamos VCF file")
    parser.add_argument("-sample", required=True, help="Name of the sample.")
    parser.add_argument("-output", required=True, help="Prefix for output files.")

    options = parser.parse_args()

    run(**vars(options))