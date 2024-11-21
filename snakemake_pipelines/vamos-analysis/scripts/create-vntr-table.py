import argparse
import sys
import numpy as np

def run(stats=None, sites=None):
    
    data = {}
    with open(sites, 'r') as sitesfile:
        for line in sitesfile:
            line = line.rstrip().split('\t')
            ru_lengths = [len(i) for i in line[3].split(',')]
            # storing all the sites as keys for the RU count and bp stats
            data[(line[0], line[1])] = [[], [], [], int(line[2]), ru_lengths]
    
    # reading individual stat files
    count=0
    for stat_file in stats.split(','):
        with open(stat_file, 'r') as file:
            count += 1
            print(f'Reading File Number {count}', file=sys.stderr)
            for line in file:
                if line.startswith('#'):
                    continue
                line = line.rstrip().split('\t')
                num_rus = [int(i) for i in line[2].split(',') if i != '.']
                bps = [int(i) for i in line[3].split(',') if i != '.']
                vntrs = [i for i in line[4].split(',') if i != '.']
                try:
                    data[(line[0], line[1])][0].extend(num_rus)
                    data[(line[0], line[1])][1].extend(bps)
                    data[(line[0], line[1])][2].extend(vntrs)
                except KeyError:
                    print(f'Did not find CHROM: {line[0]}, POS: {line[1]}', file=sys.stderr)
    
    # writing summary stats
    print('#CHROM\tREF_START\tREF_END\tNUM_HAPS\tRU_LEN_AVG\tRU_LEN_STD\tNUM_UNIQUE_VNTRS\tMAX_RUS\tMIN_RUS\tMEDIAN_RUS\t1%_RUS\t5%_RUS\t25%_RUS\t75%_RUS\t95%_RUS\t99%_RUS\tUNIQUE_NUM_RUS\tMAX_BPS\tMIN_BPS\tMEDIAN_BPS\t25%_BPS\t75%_BPS')
    for site, value in data.items():
        #print(value, file=sys.stderr)
        n_samples = len(value[0])
        if n_samples == 0:
            print(f"{site[0]}\t{site[1]}\t{value[3]}\t{n_samples}\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN")    
            continue
        ru_len = np.array(value[4])
        rus = np.array(sorted(value[0]))
        unique_rus_num = len(set(value[0]))
        unique_vntrs = len(set(value[2]))
        bps = np.array(sorted(value[1]))
        print(f"{site[0]}\t{site[1]}\t{value[3]}\t{n_samples}\t{np.mean(ru_len)}\t{np.std(ru_len)}\t{unique_vntrs}\t{int(rus.max())}\t{int(rus.min())}\t{np.median(rus)}\t{np.percentile(rus, 1)}\t{np.percentile(rus, 5)}\t{np.percentile(rus, 25)}\t{np.percentile(rus, 75)}\t{np.percentile(rus, 95)}\t{np.percentile(rus, 99)}\t{unique_rus_num}\t{int(bps.max())}\t{int(bps.min())}\t{np.median(bps)}\t{np.percentile(bps, 25)}\t{np.percentile(bps, 75)}")


if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='create-vntr-stats.py', description="Creates VNTR table from sample-wise stats file")
    parser.add_argument("-stats", required=True, help="comma-separated list of stats files")
    parser.add_argument("-sites", required=True, help="VNTR sites list")
    
    options = parser.parse_args()

    run(**vars(options))