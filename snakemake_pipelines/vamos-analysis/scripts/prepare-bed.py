import argparse
import sys

def run(ref = None, sample = None):
    
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
            line = line.strip().split('\t')
            chr = line[0]
            start = int(line[1])
            rus = None
            altanno_only_ref = False
            ru_count_only_ref = False
            altanno_h1 = None
            altanno_h2 = None
            try:
                ref_altanno, ref_end, ref_rus = ref_dict[(chr, start)]
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
            # checking if only reference vntr is present, based on the altanno sequence
            if altanno_h1 == ref_altanno:
                if (altanno_h2 is None) or ((altanno_h2 is not None) and (altanno_h2 == ref_altanno)):
                    altanno_count_ref += 1
                    altanno_only_ref = True
            # checking if only reference vntr is present, based on the ru count in the altanno
            if len(altanno_h1) == len(ref_altanno):
                if (altanno_h2 is None) or ((altanno_h2 is not None) and (len(altanno_h2) == len(ref_altanno))):
                    ru_count_ref += 1
                    ru_count_only_ref = True
            total += 1
            print(f'{chr}\t{start}\t{ref_end}\t{0 if altanno_only_ref else 1}\t{0 if ru_count_only_ref else 1}')

    print(f'#1\tNumber of sample VNTRs not found in reference:\t{num_sample_skipped}', file=sys.stderr)
    print(f'#2\tNumber of reference sites (based on ALTANNO sequence):\t{altanno_count_ref}', file=sys.stderr)
    print(f'#3\tNumber of non-reference sites (based on ALTANNO sequence):\t{total-altanno_count_ref}', file=sys.stderr)
    print(f'#4\tNumber of reference sites (based on RU Count in ALTANNO):\t{ru_count_ref}', file=sys.stderr)
    print(f'#5\tNumber of non-reference sites (based on RU Count in ALTANNO):\t{total-ru_count_ref}', file=sys.stderr)

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='prepare-bed.py', description="Processes reference VNTRs and sample VNTRs and output site list with non-reference sites.")
    parser.add_argument("-ref", required=True, help="Reference VNTRs (made from contigs)")
    parser.add_argument("-sample", required=True, help="Sample VNTRs (made from reads)")

    options = parser.parse_args()

    run(**vars(options))