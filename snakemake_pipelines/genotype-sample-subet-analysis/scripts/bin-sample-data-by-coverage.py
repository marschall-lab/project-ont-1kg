import argparse
import sys
import pandas as pd

def run(binsize = None, tsv = None, output = None):
    
    cov_data = pd.read_csv(tsv, sep='\t', header=0)
    cov_data = cov_data[['SAMPLE', 'T2T_median_cov']].rename(columns={'SAMPLE': 'sample', 'T2T_median_cov': 'cov'}).sort_values(by=['cov'])
    
    max_cov = cov_data['cov'].max()
    x=0
    while x < max_cov:
        samples = cov_data.loc[(cov_data['cov'] < x+binsize) & (cov_data['cov'] >= x)]['sample'].to_list()
        print(f'Coverages: {x}-{x+binsize-1}\tNumber of Samples: {len(samples)}', file=sys.stderr)
        with open(f'{output}_{x}-{x+binsize-1}.tsv', 'w') as outfile:
            for s in samples:
                print(s, file=outfile)
        x += binsize
    

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='bin-sample-data-by-coverage.py', description="bin samples by coverages")
    parser.add_argument("-binsize", required=True, type=int, help="Binsize of median coverage")
    parser.add_argument("-tsv", required=True, help="TSV file with the sample coverage")
    parser.add_argument("-output", required=True, help="Output prefix for sample lists")

    options = parser.parse_args()

    run(**vars(options))