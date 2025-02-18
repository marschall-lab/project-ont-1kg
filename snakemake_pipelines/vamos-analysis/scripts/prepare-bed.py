import argparse
import sys

def run(ref = None, sample = None):
    
    # parse reference VNTRs
    ref_dict = {}
    num_ref_len_mismatch = 0
    with open(ref, 'r') as vcffile:
        for line in vcffile:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            chr = line[0]
            start = int(line[1])
            end = None
            rus = None
            len_match = False
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
            altanno_len = 0
            for allele in altanno:
                altanno_len += len(rus[int(allele)])
            if altanno_len == (end-start):
                len_match = True
            else:
                num_ref_len_mismatch += 1
            ref_dict[(chr, start)] = [altanno, altanno_len, end, len_match, rus]
    
    # parse sample VNTRs
    num_sample_skipped = 0
    with open(sample, 'r') as vcffile:
        for line in vcffile:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            chr = line[0]
            start = int(line[1])
            rus = None
            only_ref = False
            altanno_h1 = None
            altanno_h2 = None
            try:
                ref_altanno, ref_altanno_len, ref_end, ref_len_match, ref_rus = ref_dict[(chr, start)]
            except KeyError:
                num_sample_skipped += 1
                continue
            info = line[7]
            for field in info.split(';'):
                if field.startswith('RU'):
                    rus = tuple(field.split('=')[1].split(','))
                if field.startswith('ALTANNO_H1'):
                    altanno_h1 = tuple(field.split('=')[1].split(','))
                if field.startswith('ALTANNO_H2'):
                    altanno_h2 = tuple(field.split('=')[1].split(','))
            assert rus == ref_rus
            # checking if only reference vntr is present
            if altanno_h1 == ref_altanno:
                if (altanno_h2 is None) or ((altanno_h2 is not None) and (altanno_h2 == ref_altanno)):
                    only_ref = True
            if only_ref:
                continue
            print(f'{chr}\t{start}\t{ref_end}\t{0 if ref_len_match else 1}\t{ref_altanno_len}')

    print(f'Number of VNTRs where reference length and vamos allele length does not match: {num_ref_len_mismatch}', file=sys.stderr)
    print(f'Number of sample VNTRs not found in reference: {num_sample_skipped}', file=sys.stderr)

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='prepare-bed.py', description="Processes reference VNTRs and sample VNTRs and output site list with non-reference sites.")
    parser.add_argument("-ref", required=True, help="Reference VNTRs (made from contigs)")
    parser.add_argument("-sample", required=True, help="Sample VNTRs (made from reads)")

    options = parser.parse_args()

    run(**vars(options))