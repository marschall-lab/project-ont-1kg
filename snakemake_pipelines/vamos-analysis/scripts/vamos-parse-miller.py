import argparse
import sys
from collections import Counter, defaultdict

def run(vcf = None, sample = None, output = None):
    
    count_writer = open(output+'-count.tsv', 'w')
    seq_writer = open(output+'-sequence.tsv', 'w')

    sample_to_idx = []
    # counting number of heterozygous VNTRs
    het_count = 0
            
    # vcf should not be gzipped
    with open(vcf, 'r') as vcffile:
        for line in vcffile:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                samples = line.rstrip().split('\t')[9:]
                #n_samples = len(list(set(samples)))
                #print(f"Number of Samples = {n_samples}")
                for idx, s in enumerate(samples):
                    if sample != s:
                        continue
                    sample_to_idx.append(9+idx)
                continue
            line = line.rstrip().split('\t')
            chrom = line[0]
            pos = line[1]
            info_field = line[7]
            rus = None
            altanno = None
            af = defaultdict(lambda: 0)     # map between genotype and the number of times the allele is found.
            tot = 0
            for info in info_field.split(';'):
                # reading info on the repeating units
                if info.startswith('RU='):
                    rus = info.split('=')[1].split(',')
                # reading info on the alternate annotation
                if info.startswith('ALTANNO='):
                    assert rus != None
                    altanno = info.split('=')[1].split(',')
            for gts in line[9:]:
                if gts[0] == '.':
                    continue
                gt = int(gts.split('/')[0])
                af[gt-1] += 1
                tot += 1
                pass
            #for index in range(n_samples):
            #    gt1 = line[9+(2*index)].split('/')[0]
            #    gt2 = line[10+(2*index)].split('/')[0]
            #    if gt1 != gt2:
            #        het_count += 1
            counter = Counter()
            for index in sample_to_idx:
                gt = line[index].split('/')[0]
                if gt == '.':
                    continue
                seq=''.join([rus[int(x)] for x in altanno[int(gt)-1].split('-')])
                print(f"{chrom}\t{pos}\t{seq}\t{gt}\t{af[int(gt)-1]/tot}", file=seq_writer)
                counter.update([int(x) for x in altanno[int(gt)-1].split('-')])
           
            for ru_idx, ru in enumerate(rus):
                if ru_idx in counter:
                    print(f"{chrom}\t{pos}\t{ru}\t{counter[ru_idx]}", file=count_writer)
                else:
                    print(f"{chrom}\t{pos}\t{ru}\t0", file=count_writer)
    
    #print(f"Number of HET VNTRs = {het_count}")
    count_writer.close()
    seq_writer.close()

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='vamos-parse-miller.py', description="Collects the repeating units from the Miller VCF and converts it into a TSV with counts")
    parser.add_argument("-vcf", required=True, help="Vamos VCF file")
    parser.add_argument("-sample", required=True, help="Name of the sample.")
    parser.add_argument("-output", required=True, help="Prefix for output files.")

    options = parser.parse_args()

    run(**vars(options))