import argparse

def convert_to_seq(anno, rus):
    anno = anno.split('-')
    seq = ''
    for a in anno:
        seq += rus[int(a)]
    return seq

def run(vcf=None, sample=None, reference=None):

    reference_dict = {}
    reader = open(reference, 'r')
    for line in reader:
        chrom, start, end, seq = line.strip().split('\t')
        reference_dict[(chrom, start, end)] = seq
    reader.close()
    
    reader = open(vcf, 'r')
    for line in reader:
        line = line.strip()
        if line.startswith('##'):
            if "<ID=SVTYPE," in line:
                continue
            if "<ID=SVLEN," in line:
                continue
            if "<ID=VNTR," in line:
                continue
            print(line)
            continue
        if line.startswith('#'):
            line = line.split('\t')
            hap1_index = line.index(sample+'_hp1')
            hap2_index = line.index(sample+'_hp2')
            header = line[0:9]+[sample]
            print('\t'.join(header))
            continue
        line = line.split('\t')
        chrom = line[0]
        start = line[1]
        info_dict = {}
        for field in line[7].split(';'):
            field = field.split('=')
            assert len(field) == 2
            info_dict[field[0]] = field[1]
        assert 'RU' in info_dict
        assert 'ALTANNO' in info_dict
        assert 'END' in info_dict
        info_field = ';'.join([f'{s}={info_dict[s]}' for s in ['END', 'RU', 'ALTANNO']])
        end = info_dict['END']
        ref_seq = reference_dict[(chrom, start, end)]
        rus = info_dict['RU'].split(',')
        altanno = info_dict['ALTANNO'].split(',')
        gt1 = line[hap1_index].split('/')[0]
        gt2 = line[hap2_index].split('/')[0]
        if gt1 == '.' and gt2 == '.':
            continue
        if gt1 == '.':
            gt2 = int(gt2)
            alt_allele_2 = convert_to_seq(altanno[gt2-1], rus)
            if alt_allele_2 == ref_seq:
                continue
            else:
                print(f'{chrom}\t{start}\t.\t{ref_seq}\t{alt_allele_2}\t.\tPASS\t{info_field}\tGT\t1/0')
                continue
        if gt2 == '.':
            gt1 = int(gt1)
            alt_allele_1 = convert_to_seq(altanno[gt1-1], rus)
            if alt_allele_1 == ref_seq:
                continue
            else:
                print(f'{chrom}\t{start}\t.\t{ref_seq}\t{alt_allele_1}\t.\tPASS\t{info_field}\tGT\t1/0')
                continue
        
        alt_allele_1 = convert_to_seq(altanno[int(gt1)-1], rus)
        alt_allele_2 = convert_to_seq(altanno[int(gt2)-1], rus)
        if (alt_allele_1 != ref_seq) and (alt_allele_2 != ref_seq):
            if (alt_allele_1 == alt_allele_2):
                print(f'{chrom}\t{start}\t.\t{ref_seq}\t{alt_allele_1}\t.\tPASS\t{info_field}\tGT\t1/1')
            else:
                # lines have to be split
                print(f'{chrom}\t{start}\t.\t{ref_seq}\t{alt_allele_1}\t.\tPASS\t{info_field}\tGT\t1/0')
                print(f'{chrom}\t{start}\t.\t{ref_seq}\t{alt_allele_2}\t.\tPASS\t{info_field}\tGT\t1/0')
            continue
        if (alt_allele_1 == ref_seq) and (alt_allele_2 == ref_seq):
            #skip line since both are reference
            continue
        if (alt_allele_1 == ref_seq):
            print(f'{chrom}\t{start}\t.\t{ref_seq}\t{alt_allele_2}\t.\tPASS\t{info_field}\tGT\t1/0')
        else:
            print(f'{chrom}\t{start}\t.\t{ref_seq}\t{alt_allele_1}\t.\tPASS\t{info_field}\tGT\t1/0')
    reader.close()
        
    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='subsample-vcf.py', description="subsample miller lab's vamos vcf and makes it compatible with the callset comparision pipeline")
    parser.add_argument("-vcf", required=True, help="Miller Vamos VCF")
    parser.add_argument("-sample", required=True, help="Sample to extract")
    parser.add_argument("-reference", required=True, help="Reference BED file")

    options = parser.parse_args()

    run(**vars(options))