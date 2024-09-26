import argparse
import sys
import pandas as pd

def run(binsize = None, tsv = None, output = None):
    
    n50_data = pd.read_csv(tsv, sep='\t', header=0)
    n50_data = n50_data[['SAMPLE', 'READ_N50']].rename(columns={'SAMPLE': 'sample', 'READ_N50': 'n50'}).sort_values(by=['n50'])
    
    max_n50 = n50_data['n50'].max()
    x=0
    while x < max_n50:
        samples = n50_data.loc[(n50_data['n50'] < x+binsize) & (n50_data['n50'] >= x)]['sample'].to_list()
        print(f'Coverages: {x}-{x+binsize-1}\tNumber of Samples: {len(samples)}', file=sys.stderr)
        with open(f'{output}_{x}-{x+binsize-1}.tsv', 'w') as outfile:
            for s in samples:
                print(s, file=outfile)
        x += binsize
    

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='bin-sample-data-by-n50.py', description="bin samples by n50s")
    parser.add_argument("-binsize", required=True, type=int, help="Binsize of n50")
    parser.add_argument("-tsv", required=True, help="TSV file with the sample n50")
    parser.add_argument("-output", required=True, help="Output prefix for sample lists")

    options = parser.parse_args()

    run(**vars(options))