import argparse
import sys
import pandas

def run(stats=None, sites=None):
    
    data = {}
    with open(sites, 'r') as sitesfile:
        for line in sitesfile:
            line = line.rstrip().split('\t')
            # storing all the sites as keys for the RU count and bp stats
            data[(line[0], line[1])] = [[], [], [], line[2]]
    
    # reading individual stat files
    for stat_file in stats.split(','):
        with open(stat_file, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                line = line.rstrip().split('\t')
                num_rus = [int(i) for i in line[2].split(',')]
                bps = [int(i) for i in line[3].split(',')]
                vntrs = line[4].split(',')
                try:
                    data[(line[0], line[1])][0].extend(num_rus)
                    data[(line[0], line[1])][1].extend(bps)
                    data[(line[0], line[1])][2].extend(vntrs)
                except KeyError:
                    print(f'Did not find CHROM: {line[0]}, POS: {line[1]}', file=sys.stderr)
    
    # writing summary stats
    print('#CHROM\tREF_START\tREF_END\tNUM_SAMPLES\tNUM_UNIQUE_VNTRS\tMAX_RUS\tMIN_RUS\tMEDIAN_RUS\t25%_RUS\t75%_RUS\tUNIQUE_NUM_RUS\tMAX_BPS\tMIN_BPS\tMEDIAN_BPS\t25%_BPS\t75%_BPS')
    for site, value in data.items():
        n_samples = len(value[0])
        if n_samples == 0:
            print(f"{site[0]}\t{site[1]}\t{value[3]}\t{n_samples}\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN")    
            continue
        rus = pandas.Series(value[0]).describe()
        unique_rus_num = len(set(value[0]))
        unique_vntrs = len(set(value[2]))
        bps = pandas.Series(value[1]).describe()
        print(f"{site[0]}\t{site[1]}\t{value[3]}\t{n_samples}\t{unique_vntrs}\t{int(rus['max'])}\t{int(rus['min'])}\t{rus['50%']}\t{rus['25%']}\t{rus['75%']}\t{unique_rus_num}\t{int(bps['max'])}\t{int(bps['min'])}\t{bps['50%']}\t{bps['25%']}\t{bps['75%']}")


if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='create-vntr-stats.py', description="Creates VNTR table from sample-wise stats file")
    parser.add_argument("-stats", required=True, help="comma-separated list of stats files")
    parser.add_argument("-sites", required=True, help="VNTR sites list")
    
    options = parser.parse_args()

    run(**vars(options))