import argparse
from collections import defaultdict
import sys

def run(bed=None, map=None, vcf=None):
    
    bed_reader = open(bed, 'r')
    aa_dict = defaultdict(lambda: None)
    for line in bed_reader:
        if line[0] == '#':
            continue
        _, _, _, bub_id, is_alt, _, allele_path = line.strip().split('\t')
        aa_dict[bub_id] = [is_alt, allele_path]
    bed_reader.close()
    
    map_reader = open(map, 'r')
    bub_to_sv_map = defaultdict(lambda: None)
    sv_to_bub_map = defaultdict(lambda: None)
    for line in map_reader:
        if line[0] == '#':
            continue
        bub_id, sv_ids = line.strip().split('\t')
        sv_ids = sv_ids.split(',')
        bub_to_sv_map[bub_id] = sv_ids
        for sv_id in sv_ids:
            sv_to_bub_map[sv_id] = bub_id
    map_reader.close()

    stat_counter = defaultdict(lambda: 0)
    vcfreader = open(vcf, 'r')
    for line in vcfreader:
        if line.startswith('##'):
            print(line.strip())
            continue
        if line.startswith('#'):
            # write the info line for the ancestral allele
            print('##INFO=<ID=ANCESTRAL_ALLELE,Number=1,Type=String,Description="Ancestral Allele Annotation. \'.\' = unknown ancestral allele | \'0\' = reference ancestral allele | \'1\' = alternate ancestral allele | \'2\' = ancestral allele is alternate but does not match the sv allele">')
            print(line.strip())
            continue
        line = line.strip().split('\t')
        sv_id = line[2]
        aa = aa_dict[sv_to_bub_map[sv_id]]
        # if the ancestral allele is unknown or the reference, then
        # update the INFO field to include that information.
        if aa[0] == '0' or aa[0] == '.':
            stat_counter[f'Number of SVs with AA is {aa[0]}'] += 1
            line[7] += f';ANCESTRAL_ALLELE={aa[0]}'
            print('\t'.join(line))
            continue
        assert aa[0] == '1'
        assert aa[1] != '.'
        sv_path = sv_id.spilt('-')[3]
        aa_path = aa[1]
        if sv_path in aa_path:
            # the bubble path of the sv allele is found in the ancestral allele path
            line[7] += f';ANCESTRAL_ALLELE={aa[0]}'
            if sv_path == aa_path:
                stat_counter[f'Number of SVs with AA is {aa[0]} (Exact match)'] += 1
            else:
                stat_counter[f'Number of SVs with AA is {aa[0]} (Partial match)'] += 1
            print('\t'.join(line))
            continue
        # only remaining case is that the ancestral allele is an alternate allele but does not have the sv bubble path
        line[7] += f';ANCESTRAL_ALLELE=2'
        stat_counter[f'Number of SVs with AA is 2'] += 1
        print('\t'.join(line))
    vcfreader.close()

    for key, value in stat_counter.items():
        print(f'{key} = {value}', file=sys.stderr)
    

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='match-ancestral-allele.py', description="match the ancestral alleles to the sv alleles and annotate the vcf.")
    parser.add_argument("-bed", required=True, help="BED file with ancestral allele")
    parser.add_argument("-map", required=True, help="map between SV alleles and the bubble they are in")
    parser.add_argument("-vcf", required=True, help="Phased Final VCF")
    
    options = parser.parse_args()

    run(**vars(options))