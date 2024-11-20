import argparse
import sys

def parse_info(line):
    info_fields={}
    for field in line.split(';')[:-1]:
        assert '=' in field
        key, value = field.split('=')
        info_fields[key] = value
    return info_fields

def run(sites=None, vcfs=None):
    
    # reading sample list
    sample_list = []
    sample_to_index = {}
    for index, vcf_file in enumerate(vcfs.split(',')):
        sample_name = vcf_file.strip().split('/')[-1][0:7]
        sample_list.append(sample_name)
        sample_to_index[sample_name] = index
    
    # reading the sites list
    common_info = {}
    altanno_info = {}
    gt_info = {}
    with open(sites, 'r') as sitesfile:
        for line in sitesfile:
            line = line.rstrip().split('\t')
            key = f'{line[0]}_{line[1]}'
            # storing all the sites as keys and END and RU as values
            common_info[key] = [line[0], line[1], line[2], line[3]]
            altanno_info[key] = []
            # initialising gt info to dots
            gt_info[key] = ['./.' for _ in sample_list]
    
    is_header_written = None
    count = 0
    # reading individual vcf files
    for vcf_file in vcfs.split(','):
        sample_name = vcf_file.strip().split('/')[-1][0:7]
        count += 1
        print(f'Reading file {count}', file=sys.stderr)
        with open(vcf_file, 'r') as file:
            for line in file:
                if line.startswith('##'):
                    # skipping INFO field lines for ALTANNO and LEN
                    if line.startswith('##INFO=<ID=ALTANNO') or line.startswith('##INFO=<ID=LEN'):
                        continue
                    if is_header_written == None:
                        # printing header
                        print(line.strip())
                    continue
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                chrom = line[0]
                pos = line[1]
                info_fields = parse_info(line[7])
                key = chrom+'_'+pos
                assert '.' not in line[9]
                # checking if the common information is consistent over all vcfs
                assert common_info[key] == [chrom, pos, info_fields['END'], info_fields['RU']]
                # extracting the altanno assigning gt numbers
                alt1, alt2, gt1, gt1 = None, None, None, None
                alt1 = '-'.join(info_fields['ALTANNO_H1'].split(','))
                # checking if altanno_h1 is already in list of altannos
                if alt1 in altanno_info[key]:
                    # if present, get the gt number
                    gt1 = altanno_info[key].index(alt1)
                else:
                    # if absent, add it to list and gt number is assigned
                    altanno_info[key].append(alt1)
                    gt1 = len(altanno_info[key])-1
                # checking if the sample has altanno_h2
                if 'ALTANNO_H2' in info_fields:
                    # if present, do the same as altanno_h1
                    alt2 = '-'.join(info_fields['ALTANNO_H2'].split(','))
                    if alt2 in altanno_info[key]:
                        gt2 = altanno_info[key].index(alt1)
                    else:
                        altanno_info[key].append(alt2)
                        gt2 = len(altanno_info[key])-1
                else:
                    # if absent, altanno_h2 is same as h1
                    gt2 = gt1
                # storing the gt information
                gt_info[key][sample_to_index[sample_name]] = f'{gt1}/{gt2}'
            is_header_written = True
    print('##INFO=<ID=ALTANNO,Number=A,Type=String,Description="Motif representation for the alternate alleles">')
    print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+"\t".join(sample_list))

    # writing the output records
    for key in common_info.keys():
        chrom, pos, end, ru = common_info[key]
        altannos = ','.join(altanno_info[key])
        gts = '\t'.join(gt_info[key])
        print(f'{chrom}\t{pos}\t.\tN\t<VNTR>\t.\tPASS\tEND={end};RU={ru};SVTYPE=VNTR;ALTANNO={altannos}\tGT\t{gts}')

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='combine-vamos-vcf.py', description="Creates multisample vcf from sample vcfs")
    parser.add_argument("-sites", required=True, help="sites list")
    parser.add_argument("-vcfs", required=True, help="comma-separated list of vcf files")
    
    options = parser.parse_args()

    run(**vars(options))