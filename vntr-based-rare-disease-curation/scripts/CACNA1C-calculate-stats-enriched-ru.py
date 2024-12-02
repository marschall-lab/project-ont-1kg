'''
Script to find the statistics of the repeat units (reported as enriched in https://doi.org/10.1016/j.ajhg.2018.07.011)
in our VNTR data.

Usage: python CACNA1C-calculate-stats-enriched-ru.py -vcf vamos-multisample.vcf

vamos-multisample.vcf created by rule vamos_t2t_vcf_combine of snakemake_pipelines/vamos-analysis/rules/vamos.smk
'''
import argparse
from scipy import stats
import numpy

def parse_info(line):
    info_fields={}
    for field in line.split(';'):
        assert '=' in field
        key, value = field.split('=')
        info_fields[key] = value
    return info_fields

def run(vcf):

    #original_motifs = 'GACCCTGACCTGACTAGTTTACAATCACAC,GATCCTGACCTGACTAGTTTACAATCACAC,GACCCTGACCTGACTAGTTTACAACCACAC,GATCCTGACCTTGCTAGTTTACAATCACAC,TATCCTGACCTGACTAGTTTACAATCACAC,GACCCTGACCTGACTAGTTTACAGTCACAC,GACCCTGACCTTACTAGTTTACAATCACAC,GATCCTGACCTGACTAGTTTACAACCACAC,GATCCTGACCTTACTAGTTTACAATCACAC,GATCCTGACCTTACTAGTTTACAATCACAG,GACCCTGACCTTACTAGTTTACAACCACAC,GATCCTGACCTTGCTAGTTTACAACCACAC,GACCCTGACCTGACTAGTTTACGATCACAC,TATCCTGACCTGACTAGTTTACAACCACAC,GACCCTGACCTGACTAGTTTACGACCACAC,GACCCTGACCTTGCTAGTTTACAATCACAC,GACCCTGACCTTACTAGTTTACGATCACAC,GATCCTGACCTTGCTAGTTTACAATCACAG,GATCCTGACCTTACTAGTTTACAATCACAA,AATCCTGACCTTGCTAGTTTACAATCACAC,GATCCTGACCTGACTACTTTACAATCACAC,TATCCTGACCTGACTAGTTTACAATCACAA,TATCCTGACCTTACTAGTTTACAATCACAA,GACCCTGACCTGACTAGTTTACAATCACAT,GATCCTGACCTGACTAGTTTACAGTCACAT,GATCCTGACCTGACTAGTTTACAGTCACAC,ACCCTGACCTTACTAGTTTACAATCACAC,GACCCTGACCTGACTAGTTTACGATCACAT,GATCCTGACCTTACTAGTTTACAATCACA,GATCCTGACCTGACTAGTTTACAATCACAT,GACCCTGACCTGACTAGTTTACAATCACAG,GACCCTGACGTGACTAGTTTACAACCACAC,GACCCTGACCTTGCTAGTTTACAACCACAC,GACCCTGACCTTACTAGTTTACGACCACAC,GATCCTGACCTGACTAGTTTACAATCACAG,GATCCTGACCTTGCTAGTTTACAATCACAA,GATCCTGACCTGACTAGTTTACGATCACAC,GACCCTGATCTGACTAGTTTACAACCACAC,GATCCTGACCTTGCTGGTTTACAACCACAC,GATCCTGACCTGACTAGTTTACGACCACAC,GATCCTGACCTTGCTGGTTTACAATCACAC,GACCCTGACCTGAGTAGTTTACAATCACAC,GACCCTGACCTGACTAGTTTACAACCAGAC,GATCCTGACGTGACTAGTTTACAACCACAC,GACCCTGACCTGGCTAGTTTACAATCACAC,GACCCCGACCTGACTAGTTTACAACCACAC,GATCCTTACCTGACTAGTTTACAATCACAC,GACCCTGACCTGACTAGTTTACGACCACAA,GACCCTGACCGGACTAGTTTACAACCACAC,GACCCTGACCTGACTAGTATACAACCACAC,GACCCTGACCTTGCTAGTTTACAACCACAA,GATCCTGACCTGACTAGTTTACCACCACAC,TTCCCTGACCTGACTAGTTTACCACCACAC,TACCCTGACTTGACTAGTTTACAATCACAC,GACCCTGACCTTGCTAGTTTACAAACACAC,GACCCTGACCTGACTAGTTTACAACCACAT,GACCCTGATCTGACTAGTTTACAATCACAC,GACCCTGACCTGACTCGTTTACAATCACAC,GACCCTGATCTTACTAGTTTACAATCACAC,GATCCTGACCTTGCTAGTTTACAACCACAA,GACCGTGACCTGACTAGTTTACAACCACAC,GACCCTGACCTGACTAGTTTACAATCACAA,GATCCTGACCTTGCTAGTTTAAAATCACAC,GACCCTGACCTGACTAGTTTAGAATCACAC,GACCCTGACCTGACTAGTTTACAAGCACAC,GACCCTGACCTGACTAGTTTACATCCACAC,GACCCTGACCTGACTTGTTTACAATCACAC,GATCCTAACCTTGCTAGTTTACAATCACAC,GACCCTCACCTGACTCGTTTACAACCACAC,GACCCTGACCTGACTAGTTTACAATAACAC,GACCCTGAACTGACTAGTTTACAACCACAC,AACCCTGACCTGACTAGTTTACAACCACAC,GACCCTGACCTGACTGGTTTACAACCACAC,GACCCTGACCTGACTAGTTTACAATTACAC,GTCCCTGACCTGACTAGTTTACAACCACAC,GACCCTGACTTGACTAGTTTACAATCACAC,GACCTTGACCTGACTAGTTTACAATCACAC,GACGCTGACCTGACTAGTTTACAATCACAC,GACCCTGACCTGACTAGTTTACAACCCGC,GACCCTGACCTGACTAGTTTACAGTCACAT,GATCCTGACCTTGGTAGTTTACAATCACAC,TATCCCGACCTGACTAGTTTACAATCACAA,AACCCTGACCTGACTAGTTTACAATCACAC,GACCCTGACCTGACTAGTTCACAACCACAC,GACCCTGACCTGACTAGTTTACAACCACAA,GACCCTGACCTGACTAGTTTACAACCACAG,GACCCTGACCTGACTAGTTTACAATCTCAC,GACGCTGACCTGACTAGTTTACAACCACAC,GACCCTGACCTGACTAGGTTACAACCACAC,GACACTGACCTGACTAGTTTACAATCACAC,GACCCTGAGCTGACTAGTTTACAATCACAC,GAACCTGACCTTGCTAGTTTACAATCACAC,GACCGTGACCTGACTAGTTTACAATCACAC,GACTCTGACCTGACTAGTTTACAATCACAC,CATCCTGACCTGACTAGTTTACAATCACAC,GACCCTGACCTGACTAGTTTACATTCACAC,GACCCTGACCTGACTAGTTTACAATCACCC,TATCCTGACCTGACTAGTTTACAATCACAT,GACCCTGACCTTACTAGTTTACAATCACAG,GATCCTGACCTTACTAGGTTACAATCACAC,GACCCTGACCTGACTAGTTTAAAATCACAC,GATCCTGACCTTGATAGTTTACAACCACAC,GATCCTTACCTGACTAGTTTACGACCACAC,GATCCTGACCTGACTAGTTTACAACCACAT,GATCCTGACTTTACTAGTTTACAATCACAG,GATCCTGACCTGACCAGTTTACAATCACAC,GACCCGGACCTGACTAGTTTACAACCACAC,GACCCTGACCTGAATAGTTTACAACCACAC,GACACTGACCTGACTAGTTTACAACCACAC,GACCCTGACCTGACGAGTTTACAATCACAC,GACCCTGACGTGACTAGTTTACAATCACAC,GACCCTGACCTGACTACTTTACAATCACAC,GATCCTGACCTGACTAGTTTACAATCACAA,GACCCTGACCTGACTAGTTTACACTCACAC,GACCCTGACCTGACTAGTTTACAACCGCAC,GATCCTGACCTGACTAGTTTACAACCACAG,GACCCTGACCTGACTAGTTTACGATCACAG,GACCCTGACCTTACTAGTTTACGATCACAT,GATCCTGACCTTACTAGTTTACAGTCACAA,AATCCTGACCTGACTAGTTTACAATCACAC,GACCCTGACCTTACTAGTTTACAATAACAC,GATCCTGACGTTGCTAGTTTACAATCACAC,AATCCTGACCTTACTAGGTTACAATCACAC,GATCCTGACCTGACTAGTTTAAAATCACAG,GACCCTGACCTGACTAGTTCACAATCACAC,GATCTTGACCTTGCTAGTTTACAACCACAC,TACCCTGACCTGACTAGTTTACAATCACAC,GACCCTAACCTGACTAGTTTACAATCACAC,GACCCTGACCTGACTAGTTTACAATCAGAC,GACCTTGACCTGACTAGTTTACAACCACAC,GATCCTCACCTTGCTAGTTTACAATCACAC,GACCCTGACCTGACTAGTTTACAACCTCAC,TATCCTGACCTGACTAGTTTACGATCACAC,GATCCTGACCGTGCTAGTTTACAACCACAC,GATCCTGGCCTGACTAGTTTACAATCACAC,GATCCTGACTTTGCTAGTTTACAATCACAC,GACCCTGACCTGACTAGTTTACCATCACAC,GACCCTGACCTGACTAGTTTCCAATCACAC,GATCCTGACCTTAGTAGTTTACAATTACAC,GATCCTGACCTGACTAGTTTACGATCATAC,GATCCTGACCTGACTAGTTTACAGCCACAC,GACCCTGACCTGTCTAGTTTACAATCACCC,GACCCTGATCTGACTAGTTTACAGTCACAC,TATCCTGACCTGACTAGTTTACAGTCACAC'.split(',')

    # hardcoded location and interested repeat units. Interested repeats from paper.
    chrom='chr12'
    start='2251595'

    prot_rus=['GATCCTGACCTGACTAGTTTACAATCACAC',
            'GATCCTGACCTTACTAGTTTACAATCACAG',
            'GACCCTGACCTGACTAGTTTACGACCACAC',
            'GACCCTGACCTGACTAGTTTACGATCACAT',
            'GACCCTGACCTTACTAGTTTACAACCACAC']
    risk_rus=['GACCCTGACCTTACTAGTTTACGATCACAC',
            'GACCCTGACCTGACTAGTTTACGATCACAC',
            'GACCCTGACCTTACTAGTTTACGACCACAC',
            'GATCCTGACCTGACTAGTTTACAATCACAT',
            'GACCCTGACCTGACTAGTTTACAACCACAC']
    
    prot_frac = []
    risk_frac = []
    cand_ru_frac = {}
    vcfreader = open(vcf, 'r')
    for line in vcfreader:
        if line[0] == '#':
            continue
        line = line.strip().split('\t')
        if not (line[0] == chrom and line[1] == start):
            continue
        info_fields = parse_info(line[7])
        rus = info_fields['RU'].split(',')
        #rus = original_motifs
        cand_rus_to_numbers = {}
        prot_rus_to_index = {}
        for ru in prot_rus:
            try:
                index = rus.index(ru)
                cand_rus_to_numbers[ru] = []
                cand_ru_frac[ru] = []
                prot_rus_to_index[ru] = index
            except ValueError:
                print(f'Protective RU {ru} not in the vcf')
        risk_rus_to_index = {}
        for ru in risk_rus:
            try:
                index = rus.index(ru)
                cand_rus_to_numbers[ru] = []
                cand_ru_frac[ru] = []
                risk_rus_to_index[ru] = index
            except ValueError:
                print(f'Risk RU {ru} not in the vcf')
        risk_index = list(risk_rus_to_index.values())
        prot_index = list(prot_rus_to_index.values())
        #print(prot_rus_to_index)
        #print(prot_index)
        #print(risk_rus_to_index)
        #print(risk_index)
        num_risk_altanno = []   # number of risk rus in altanno
        num_prot_altanno = []   # number of protective rus in altanno
        num_othr_altanno = []   # number of other rus in altanno
        for altanno in info_fields['ALTANNO'].split(','):
            rus=altanno.split('-')
            num_risk=0
            num_prot=0
            num_othr=0
            for key in cand_rus_to_numbers.keys():
                cand_rus_to_numbers[key].append(0)
            for ru in rus:
                if int(ru) in prot_index:
                    num_prot += 1
                    cand_rus_to_numbers[info_fields['RU'].split(',')[int(ru)]][-1] += 1
                elif int(ru) in risk_index:
                    num_risk += 1
                    cand_rus_to_numbers[info_fields['RU'].split(',')[int(ru)]][-1] += 1
                else:
                    num_othr += 1    
            num_prot_altanno.append(num_prot)
            num_risk_altanno.append(num_risk)
            num_othr_altanno.append(num_othr)
        #print(num_othr_altanno)
        #print(num_prot_altanno)
        #print(num_risk_altanno)
        gts=line[9:]
        for gt in gts:
            for gt in gt.split('/'):
                if gt == '.':
                    continue
                gt = int(gt)-1
                tot = num_prot_altanno[gt]+num_risk_altanno[gt]+num_othr_altanno[gt]
                prot_frac.append(num_prot_altanno[gt]/tot)
                risk_frac.append(num_risk_altanno[gt]/tot)
                for ru in cand_rus_to_numbers.keys():
                    cand_ru_frac[ru].append(cand_rus_to_numbers[ru][gt]/tot)
        break
    vcfreader.close()

    print()
    print('Stats for Protective RUs')
    print('Mean: ', numpy.mean(prot_frac))
    print('Stddev: ', numpy.std(prot_frac))
    print('\nStats for Risk RUs')
    print('Mean: ', numpy.mean(risk_frac))
    print('Stddev: ', numpy.std(risk_frac))
    
    for ru in cand_ru_frac.keys():
        title=None
        if ru in prot_rus:
            title='Protective'
        else:
            title='Risk'
        print(f'\n{ru} ({title})')
        print('Mean: ', numpy.mean(cand_ru_frac[ru]))
        print('Stddev: ', numpy.std(cand_ru_frac[ru]))
        print()

    pass

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='CACNA1C-calculate-stats-enriched-ru.py', description="Reports statistics about the RUs associated with disease")
    parser.add_argument("-vcf", required=True, help="Multisample vamos vcf")
    
    options = parser.parse_args()

    run(**vars(options))