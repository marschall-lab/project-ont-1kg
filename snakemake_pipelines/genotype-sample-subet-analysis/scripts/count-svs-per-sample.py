import argparse
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

def run(sample_data=None, callset=None, sniffles=None, delly=None, svarp=None):
    
    sample_to_pop_map = map_sample_to_pop(sample_data)
    sv_counter = {'callset': None,
                  'sniffles': None,
                  'delly': None,
                  'svarp': None}
    file_mapper = {'callset': callset,
                  'sniffles': sniffles,
                  'delly': delly,
                  'svarp': svarp}
    for vcf_type in sv_counter.keys():
        reader=open(file_mapper[vcf_type], 'r')
        sample_list = None
        for line in reader:
            if line.startswith('##'):
                continue
            if line[0] == '#':
                # count samples
                sample_list = line.strip().split('\t')[9:]
                n_samples = len(line.strip().split('\t')[9:])
                sv_counter[vcf_type] = [0 for i in range(n_samples)]
                continue
            gts = read_gts(line)
            assert len(gts) == n_samples
            sv_counter[vcf_type] = list(map(add, sv_counter[vcf_type], gts))
        reader.close()
        print(f'{vcf_type}_AFR\t{",".join([str(i) for i in sv_counter[vcf_type] if sample_to_pop_map[sample_list[i]] == 'AFR'])}')
        print(f'{vcf_type}_non_AFR\t{",".join([str(i) for i in sv_counter[vcf_type] if sample_to_pop_map[sample_list[i]] != 'AFR'])}')


    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='count-svs-per-sample.py', description="count svs per sample for some coverage run")
    parser.add_argument("-sample-data", required=True, help="Sample ancestry data")
    parser.add_argument("-callset", required=True, help="Phased callset vcf")
    parser.add_argument("-sniffles", required=True, help="Sniffles vcf")
    parser.add_argument("-delly", required=True, help="Delly vcf")
    parser.add_argument("-svarp", required=True, help="Svarp vcf")
    
    options = parser.parse_args()

    run(**vars(options))