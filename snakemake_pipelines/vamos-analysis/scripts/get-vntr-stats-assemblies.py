import argparse
from collections import defaultdict

def read_vcf(vcffile, data, index):
    for line in vcffile:
        if line.startswith('#'):
            continue
        line = line.rstrip().split('\t')
        info_field = line[7]
        chrom = line[0]
        pos = line[1]
        key = chrom+'_'+pos
        rus = None
        altanno = []
        for info in info_field.split(';'):
            if info.startswith('RU='):
                rus = info.split('=')[1].split(',')
            if info.startswith('ALTANNO_H1='):
                assert rus != None
                altanno.append(info.split('=')[1])
        num_rus = []
        bps = []
        alt_out = []
        for alt in altanno:
            alt = [int(i) for i in alt.split(',') if i != 'DEL']
            num_rus.append(str(len(alt)))
            bps.append(str(sum([len(rus[i]) for i in alt])))
            alt_out.append('-'.join([str(i) for i in alt]))
        data[key][0][index] = num_rus[0]
        data[key][1][index] = bps[0]
        data[key][2][index] = alt_out[0]

def run(hap1 = None, hap2 = None):
    
    print('#CHROM\tPOS\tNUM_RUS_VNTR\tBPS_VNTR\tVNTR')
    # to store data for Num RUs, Bps, and VNTR annots
    data = defaultdict(lambda: [[None, None] for _ in range(3)])
    with open(hap1, 'r') as vcffile:
        read_vcf(vcffile, data, 0)
    with open(hap2, 'r') as vcffile:
        read_vcf(vcffile, data, 1)
    
    for key in data.keys():
        chrom, pos = key.split('_')
        num_rus, bps, alt_out = data[key]
        for i in [0,1]:
            if num_rus[i] == None:
                num_rus[i] == ''
            if bps[i] == None:
                bps[i] == ''
            if alt_out[i] == None:
                alt_out[i] == ''
        print(f'{chrom}\t{pos}\t{",".join(num_rus)}\t{",".join(bps)}\t{",".join(alt_out)}')

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='get-vntr-stats-assemblies.py', description="Gets the VNTR stats of VAMOS calls from haplotype resolved assemblies.")
    parser.add_argument("-hap1", required=True, help="Vamos VCF file for haplotype 1")
    parser.add_argument("-hap2", required=True, help="Vamos VCF file for haplotype 2")
    
    options = parser.parse_args()

    run(**vars(options))