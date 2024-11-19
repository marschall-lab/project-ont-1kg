import argparse
import sys
from operator import add

def map_sample_to_pop(file):
    sample_to_pop = {}
    with open(file, 'r') as reader:
        for line in reader:
            if line[0] == 's':
                continue
            line = line.rstrip().split('\t')
            sample_to_pop[line[0]] = line[3]
    return sample_to_pop

def read_gts(line):
    fields = line.split('\t')[8]
    assert 'GT' in fields
    gt_index = fields.split(':').index('GT')
    gts = [i.split(':')[gt_index] for i in line.strip().split('\t')[9:]]
    counter = [0 for _ in range(len(gts))]
    for index, gt in enumerate(gts):
        if '/' in gt:
            #unphased
            separator = '/'
        elif '|' in gt:
            #phased
            separator = '|'
        else:
            # some issue with the record
            continue
        count = 0
        gt = gt.split(separator)
        assert len(gt) == 2
        for g in gt:
            if g != '.' and g != '0':
                count += 1
        counter[index] += count
    
    return counter

def run(sample_data=None, sample_list=None, path=None):
    
    sample_to_pop_map = map_sample_to_pop(sample_data)
    with open(sample_list, 'r') as sample_list_reader:
        sv_counter_list = {'AFR': [], 'non-AFR': []}
        for sample in sample_list_reader:
            sv_counter = None
            try:
                reader=open(path+sample.strip()+'.vcf', 'r')
            except FileNotFoundError:
                continue
            for line in reader:
                if line.startswith('##'):
                    continue
                if line[0] == '#':
                    # count samples
                    n_samples = len(line.strip().split('\t')[9:])
                    assert n_samples == 1
                    sv_counter = [0 for i in range(n_samples)]
                    continue
                gts = read_gts(line)
                assert len(gts) == n_samples
                sv_counter = list(map(add, sv_counter, gts))
            try:
                if sample_to_pop_map[sample.strip()] == 'AFR':
                    sv_counter_list['AFR'].append(str(sv_counter[0]))
                else:
                    sv_counter_list['non-AFR'].append(str(sv_counter[0]))
            except KeyError:
                continue
            reader.close()
        print(f'Number of AFR in single sample Delly = {len(sv_counter_list['AFR'])}', file=sys.stderr)
        print(f'Number of non AFR in single sample Delly = {len(sv_counter_list['non-AFR'])}', file=sys.stderr)
        print('single_sample_delly_AFR\t'+','.join(sv_counter_list['AFR']))
        print('single_sample_delly_non_AFR\t'+','.join(sv_counter_list['non-AFR']))


if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='count-svs-single-sample-delly.py', description="count svs from single sample delly bcfs")
    parser.add_argument("-sample-data", required=True, help="Sample ancestry data")
    parser.add_argument("-sample-list", required=True, help="List of single sample delly variants from a range")
    parser.add_argument("-path", required=True, help="Path to the VCF files")
    
    options = parser.parse_args()

    run(**vars(options))