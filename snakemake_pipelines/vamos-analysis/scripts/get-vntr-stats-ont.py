import argparse
import sys

def run(vcf = None):
    
    print('#CHROM\tPOS\tNUM_RUS_VNTR\tBPS_VNTR\tVNTR')
    with open(vcf, 'r') as vcffile:
        for line in vcffile:
            if line.startswith('#'):
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
            if len(altanno) == 1:
                altanno.append(altanno[0])
            num_rus = []
            bps = []
            alt_out = []
            for alt in altanno:
                #print(line, file=sys.stderr)
                #print(alt, file=sys.stderr)
                alt = [int(i) for i in alt.split(',') if i != 'DEL']
                num_rus.append(str(len(alt)))
                bps.append(str(sum([len(rus[i]) for i in alt])))
                alt_out.append('-'.join([str(i) for i in alt]))
            print(f'{chrom}\t{pos}\t{",".join(num_rus)}\t{",".join(bps)}\t{",".join(alt_out)}')

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='get-vntr-stats-ont.py', description="Gets the VNTR stats from ONT vcf for making the table")
    parser.add_argument("-vcf", required=True, help="Vamos VCF file")
    
    options = parser.parse_args()

    run(**vars(options))