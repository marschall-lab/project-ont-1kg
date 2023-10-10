import seaborn
import matplotlib.pyplot as plt
import pandas
import argparse


class Table:
  
    def __init__(self, stats) -> None:
        self.stats = stats
        pass
    
    def subset_panel(self, var_type=None, is_SV=False, is_not_SV=False):
        assert not (is_SV == True and is_not_SV==True) 
        if var_type:
            assert var_type in ['SNV', 'INS', 'DEL', 'COMPLEX']
        
        if var_type:
            subtable = self.stats[self.stats['varaint type'] == var_type]
        else:
            subtable = self.stats
        if is_SV:
            subtable=subtable[subtable['variant length'] >= 50]
        if is_not_SV:
            subtable=subtable[subtable['variant length'] < 50]
        columns = ['variant_id',
                'variant type',
                'variant length',
                'panel_allele_freq',
                'panel_alternative_alleles',
                'panel_total_alleles',
                'panel_unknown_alleles']
        
        return subtable[columns]
        

    def subset(self, var_type=None, pop_code='all', is_SV=False, is_not_SV=False):
        assert not (is_SV == True and is_not_SV==True) 
        if var_type:
            assert var_type in ['SNV', 'INS', 'DEL', 'COMPLEX']
        if pop_code != 'all':
            assert pop_code in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
        
        if var_type:
            subtable = self.stats[self.stats['varaint type'] == var_type]
        else:
            subtable = self.stats
        if is_SV:
            subtable=subtable[subtable['variant length'] >= 50]
        if is_not_SV:
            subtable=subtable[subtable['variant length'] < 50]
        
        columns = ['variant_id', 'variant type', 'variant length']
        columns.append('%s_allele_freq'%(pop_code))
        columns.append('%s_alternative_alleles'%(pop_code))
        columns.append('%s_total_alleles'%(pop_code))
        columns.append('%s_unknown_alleles'%(pop_code))
        columns.append('%s_heterozygosity'%(pop_code))
        columns.append('%s_heterozygous_genotypes'%(pop_code))
        columns.append('%s_homozygous_reference_genotypes'%(pop_code))
        columns.append('%s_homozygous_alternate_genotypes'%(pop_code))
        columns.append('%s_total_genotypes'%(pop_code))
        
        for q in [0,50,100,200]:
            columns.append('%s_GQ>=%s'%(pop_code,str(q)))
        
        return subtable[columns]
            



def plot_overall_hwe():
    '''
    Plots HWE curve for all the variant IDs for the entire cohort
    '''
    pass


parser = argparse.ArgumentParser(prog='plot-vcf-stats.py', description="Making plots out of the VCF callset stats.")
parser.add_argument('-table', metavar='TABLE', help='statistics table')
parser.add_argument('-output', metavar='OUTPUT', help='output directory')
args = parser.parse_args()

stats = Table(pandas.read_csv(args.table, sep='\t', header=0))