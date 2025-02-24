import argparse
import sys

def write_header(sample):
    print('##fileformat=VCFv4.2')
    print('##source=prepare-vntr-vcf.py')
    print('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">')
    print('##INFO=<ID=RU,Number=.,Type=String,Description="Repeat unit sequence">')
    print('##INFO=<ID=ALTANNO_REF,Number=1,Type=String,Description="Alternate annotation for the Reference VNTR">')
    print('##INFO=<ID=ALTANNO_H1,Number=1,Type=String,Description="Alternate annotation for the H1 VNTR">')
    print('##INFO=<ID=ALTANNO_H2,Number=1,Type=String,Description="Alternate annotation for the H2 VNTR">')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample)

def info_field_to_string(info_field):
    return ';'.join([f'{key}={value}' for key, value in info_field.items()])

def run(ref = None, sample = None):
    
    sample_name=sample.split('/')[-1].split('.')[0]
    write_header(sample_name)
    
    # parse reference VNTRs
    ref_dict = {}
    with open(ref, 'r') as vcffile:
        for line in vcffile:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            chr = line[0]
            start = int(line[1])
            end = None
            rus = None
            altanno = None
            info = line[7]
            assert line[9] == '1/1'
            for field in info.split(';'):
                if field.startswith('END'):
                    end = int(field.split('=')[1])
                if field.startswith('RU'):
                    rus = tuple(field.split('=')[1].split(','))
                if field.startswith('ALTANNO_H1'):
                    altanno = tuple(field.split('=')[1].split(','))
            ref_dict[(chr, start)] = [altanno, end, rus]
    
    # parse sample VNTRs
    num_sample_skipped = 0
    altanno_count_ref = 0
    ru_count_ref = 0
    total= 0
    with open(sample, 'r') as vcffile:
        for line in vcffile:
            if line.startswith('#'):
                continue
            info_field = {}
            outline = [None] * 10   
            line = line.strip().split('\t')
            chr = line[0]
            outline[0] = chr
            start = int(line[1])
            outline[1] = str(start)
            outline[2] = '.'        # ID
            outline[5] = '.'        # QUAL
            outline[6] = 'PASS'     # FILTER
            outline[8] = 'GT'       # FORMAT
            outline[3] = ''         # REF initializing
            outline[4] = ''         # ALT initializing
            rus = None
            sample_gt = line[9]
            assert sample_gt in ['1/1', '1/2']
            altanno_only_ref = False
            altanno_h1 = None
            altanno_h2 = None
            try:
                ref_altanno, ref_end, ref_rus = ref_dict[(chr, start)]
            except KeyError:
                num_sample_skipped += 1
                continue
            info_field['END'] = str(ref_end)
            info = line[7]
            for field in info.split(';'):
                if field.startswith('RU'):
                    rus = tuple(field.split('=')[1].split(','))
                if field.startswith('ALTANNO_H1'):
                    altanno_h1 = tuple(field.split('=')[1].split(','))
                if field.startswith('ALTANNO_H2'):
                    altanno_h2 = tuple(field.split('=')[1].split(','))
            assert rus == ref_rus
            info_field['RU'] = ','.join(rus)
            for alt in ref_altanno:
                outline[3] += rus[int(alt)]
            info_field['ALTANNO_REF'] = ','.join(ref_altanno)
            
            # checking if only reference vntr is present, based on the altanno sequence
            if altanno_h1 == ref_altanno:
                if (altanno_h2 is None) or ((altanno_h2 is not None) and (altanno_h2 == ref_altanno)):
                    altanno_count_ref += 1
                    altanno_only_ref = True
            if altanno_only_ref:
                continue
            
            if altanno_h2 is not None:
                if altanno_h1 != ref_altanno and altanno_h2 != ref_altanno:
                    info_field['ALTANNO_H1'] = ','.join(altanno_h1)
                    info_field['ALTANNO_H2'] = ','.join(altanno_h2)
                    h1_seq = ''
                    h2_seq = ''
                    for alt in altanno_h1:
                        h1_seq += rus[int(alt)]
                    for alt in altanno_h2:
                        h2_seq += rus[int(alt)]
                    if outline[3] != h1_seq and outline[3] != h2_seq:
                        outline[4] = ','.join([h1_seq, h2_seq])
                        outline[9] = '1/2'
                    elif outline[3] == h1_seq and outline[3] == h2_seq:
                        continue
                    elif outline[3] == h1_seq:
                        outline[4] = h2_seq
                        outline[9] = '1/0'
                    elif outline[3] == h2_seq:
                        outline[4] = h1_seq
                        outline[9] = '1/0'
                elif altanno_h1 != ref_altanno:
                    info_field['ALTANNO_H1'] = ','.join(altanno_h1)
                    h1_seq = ''
                    for alt in altanno_h1:
                        h1_seq += rus[int(alt)]
                    if outline[3] != h1_seq:
                        outline[4] = h1_seq
                        outline[9] = '1/0'
                    else:
                        continue
                elif altanno_h2 != ref_altanno:
                    info_field['ALTANNO_H1'] = ','.join(altanno_h2)
                    h1_seq = ''
                    for alt in altanno_h2:
                        h1_seq += rus[int(alt)]
                    if outline[3] != h1_seq:
                        outline[4] = h1_seq
                        outline[9] = '1/0'
                    else:
                        continue
            else:
                assert altanno_h1 != ref_altanno
                info_field['ALTANNO_H1'] = ','.join(altanno_h1)
                h1_seq = ''
                for alt in altanno_h1:
                    h1_seq += rus[int(alt)]
                if outline[3] != h1_seq:
                    outline[4] = h1_seq
                    outline[9] = '1/1'
                else:
                    continue
            outline[7] = info_field_to_string(info_field)
            print('\t'.join(outline))

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='prepare-vntr-vcf.py', description="Processes reference VNTRs and sample VNTRs and output site a VCF with non-reference records.")
    parser.add_argument("-ref", required=True, help="Reference VNTRs (made from contigs)")
    parser.add_argument("-sample", required=True, help="Sample VNTRs (made from reads)")

    options = parser.parse_args()

    run(**vars(options))