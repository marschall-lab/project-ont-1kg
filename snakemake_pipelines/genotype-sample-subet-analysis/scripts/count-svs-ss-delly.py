import argparse
from operator import add

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

def run(sample_list=None, path=None):
       
    with open(sample_list, 'r') as sample_list_reader:
        sv_counter_list = []
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
            sv_counter_list.append(str(sv_counter[0]))
            reader.close()
        print(','.join(sv_counter_list))


if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='count-svs-ss-delly.py', description="count svs from single sample delly bcfs")
    parser.add_argument("-sample-list", required=True, help="List of single sample delly variants")
    parser.add_argument("-path", required=True, help="Path to the VCF files")
    
    options = parser.parse_args()

    run(**vars(options))