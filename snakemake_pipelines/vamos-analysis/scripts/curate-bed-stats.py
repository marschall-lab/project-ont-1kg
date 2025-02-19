import argparse
import pandas as pd

def run(stats = None):
    
    skipped_sites = []
    count_ref_sites = []
    count_nref_sites = []
    for stat in stats.split(','):
        count_stats = {}
        with open(stat, 'r') as file:
            for line in file:
                line=line.strip().split('\t')
                count_stats[line[0]] = int(line[2])
        skipped_sites.append(count_stats['#1'])
        count_ref_sites.append(count_stats['#2']/(count_stats['#2']+count_stats['#3']))
        count_nref_sites.append(count_stats['#3']/(count_stats['#2']+count_stats['#3']))
    
    print('Skipped Sites Statistics:')
    print(pd.DataFrame(skipped_sites).describe())
    print('\nReference Sites Statistics:')
    print(pd.DataFrame(count_ref_sites).describe())
    print('\nNon-reference Sites Statistics:')
    print(pd.DataFrame(count_nref_sites).describe())
            
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='curate-bed-stats.py', description="Some basic stats from all sample VNTRs combined")
    parser.add_argument("-stats", required=True, help="comma-separated list of stats files")
    
    options = parser.parse_args()

    run(**vars(options))