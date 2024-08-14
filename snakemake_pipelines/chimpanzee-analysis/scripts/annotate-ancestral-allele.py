import argparse
from collections import defaultdict
import sys

def run(vcf = None):
    
    vcfreader = open(vcf, 'r')
    print("#chr\tref_start\tref_end\tbubble_id\is_alternate\ttancestral_allele\tancestral_allele_bubble_path")
    counts = defaultdict(lambda: 0)
    for line in vcfreader:
        if line.startswith('#'):
            continue
        chr, pos, bub_id, ref_seq, alt_seq, _, _, info, _, gt = line.rstrip().split('\t')
        at = None
        ancestral_allele = None
        ancestral_allele_path = None
        info = info.split(';')
        for field in info:
            if field.startswith('AT='):
                at = field[3:].split(',')
        assert at != None
        if gt == '1':
            ancestral_allele = alt_seq
            ancestral_allele_path = at[1]
            assert alt_seq != '.'
            counts['Number of Bubble with Alternate Ancestral Allele'] += 1
        elif gt == '0':
            ancestral_allele = ref_seq
            ancestral_allele_path = at[0]
            counts['Number of Bubble with Reference Ancestral Allele'] += 1
        else:
            assert gt == '.'
            ancestral_allele = '.'
            ancestral_allele_path = '.'
            counts['Number of Bubble with Unknown Ancestral Allele'] += 1
        print(f"{chr}\t{pos}\t{int(pos)+len(ref_seq)}\t{bub_id}\t{gt}\t{ancestral_allele}\t{ancestral_allele_path}")
    for key, value in counts.items():
        print(f"{key}: {value}", file=sys.stderr)

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='annotate-ancestral-allele.py', description="create TSV files with annotatations for the ancestral allele.")
    parser.add_argument("-vcf", required=True, help="VCF with haploid entry")
    
    options = parser.parse_args()

    run(**vars(options))