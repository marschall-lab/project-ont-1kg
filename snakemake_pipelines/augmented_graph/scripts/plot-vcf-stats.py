import seaborn
import matplotlib.pyplot as plt
import pandas
import argparse
import numpy as np
from scipy.stats import pearsonr

class Table:
  
    def __init__(self, stats) -> None:
        self._stats = stats
        pass
        
    def subset(self, var_type=None, is_SV=False, is_not_SV=False):
        assert not (is_SV == True and is_not_SV==True) 
        if var_type:
            assert var_type in ['SNV', 'INS', 'DEL', 'COMPLEX']
        if var_type:
            subtable = self._stats[self._stats['variant_type'] == var_type]
        else:
            subtable = self._stats
        if is_SV:
            subtable=subtable[subtable['variant_length'] >= 50]
        if is_not_SV:
            subtable=subtable[subtable['variant_length'] < 50]
        
        return subtable


def plot_all_hwe(data, out):
    '''
    Plots HWE curve for all the variant IDs for the entire cohort
    '''
    x_name="all_allele_freq"
    y_name="all_heterozygosity"
    title=out.split('/')[-1][0:-4]
    #g = seaborn.jointplot(data=data, x=x_name, y=y_name, kind='hex', height=10, ratio=5, ylim=(0, 1), cbar=True)
    g = seaborn.JointGrid(data=data, x=x_name, y=y_name, height=10, ratio=5, ylim=(0, 1))
    g.plot_marginals(seaborn.histplot, bins=50)
    g.plot_joint(plt.hexbin, bins='log', gridsize=50, cmap='Blues')
    g.set_axis_labels(x_name, y_name)
    x = np.linspace(0,1,10000)
    y = 2*x*(1-x)
    plt.plot(x,y, color='black', linewidth=2)
    plt.subplots_adjust(left=0, right=0.8, top=1, bottom=0)  # shrink fig so cbar is visible
    # make new ax object for the cbar
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    cbar_ax = g.fig.add_axes([.85, .1, .05, .8])  # x, y, width, height
    plt.colorbar(cax=cbar_ax)
    plt.yticks(fontsize=15)
    g.fig.suptitle(title)
    g.savefig(out)
    plt.close()
    pass

def plot_pop_hwe(data, out, pop):
    '''
    Plots HWE curve for all the variant IDs for the entire cohort
    '''
    x_name=pop+"_allele_freq"
    y_name=pop+"_heterozygosity"
    title=out.split('/')[-1][0:-4]
    g = seaborn.JointGrid(data=data, x=x_name, y=y_name, height=10, ratio=5, ylim=(0, 1))
    g.plot_marginals(seaborn.histplot, bins=50)
    g.plot_joint(plt.hexbin, bins='log', gridsize=50, cmap='Blues')
    g.set_axis_labels(x_name, y_name)
    x = np.linspace(0,1,10000)
    y = 2*x*(1-x)
    plt.plot(x,y, color='black', linewidth=2)
    plt.subplots_adjust(left=0, right=0.8, top=1, bottom=0)  # shrink fig so cbar is visible
    # make new ax object for the cbar
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    cbar_ax = g.fig.add_axes([.85, .1, .05, .8])  # x, y, width, height
    plt.colorbar(cax=cbar_ax)
    plt.yticks(fontsize=15)
    g.fig.suptitle(title)
    g.savefig(out)
    plt.close()
    pass


def plot_all_panelAF_vs_callsetAF(data, out):
    '''
    Plots panel allele freq vs callset allele freq for all the variant IDs
    '''
    x_name="panel_allele_freq"
    y_name="all_allele_freq"
    title=out.split('/')[-1][0:-4]
    print(title)
    print(pearsonr(data[x_name], data[y_name]))
    g = seaborn.JointGrid(data=data, x=x_name, y=y_name, height=10, ratio=5)
    g.plot_marginals(seaborn.histplot, bins=50)
    g.plot_joint(plt.hexbin, bins='log', gridsize=50, cmap='Blues')
    g.set_axis_labels(x_name, y_name)
    plt.subplots_adjust(left=0, right=0.8, top=1, bottom=0)  # shrink fig so cbar is visible
    # make new ax object for the cbar
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    cbar_ax = g.fig.add_axes([.85, .1, .05, .8])  # x, y, width, height
    plt.colorbar(cax=cbar_ax)
    plt.yticks(fontsize=15)
    g.fig.suptitle(title)
    g.savefig(out)
    plt.close()
    pass

def plot_pop_panelAF_vs_callsetAF(data, out, pop):
    '''
    Plots panel allele freq vs callset allele freq for all the variant IDs
    '''
    x_name="panel_allele_freq"
    y_name=pop+"_allele_freq"
    title=out.split('/')[-1][0:-4]
    print(title)
    print(pearsonr(data[x_name], data[y_name]))
    g = seaborn.JointGrid(data=data, x=x_name, y=y_name, height=10, ratio=5)
    g.plot_marginals(seaborn.histplot, bins=50)
    g.plot_joint(plt.hexbin, bins='log', gridsize=50, cmap='Blues')
    g.set_axis_labels(x_name, y_name)
    plt.subplots_adjust(left=0, right=0.8, top=1, bottom=0)  # shrink fig so cbar is visible
    # make new ax object for the cbar
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    cbar_ax = g.fig.add_axes([.85, .1, .05, .8])  # x, y, width, height
    plt.colorbar(cax=cbar_ax)
    plt.yticks(fontsize=15)
    g.fig.suptitle(title)
    g.savefig(out)
    plt.close()
    pass

def print_basic_stats(table):
    SVtable = table[table['variant_length'] >= 50]
    print("\tNumber of variants: ", table.shape[0])
    print("\tNumber of SVs: ", table[table['variant_length'] >= 50].shape[0])
    print("\tNumber of INS: ", table[table['variant_type'] == 'INS'].shape[0])
    print("\tNumber of DEL: ", table[table['variant_type'] == 'DEL'].shape[0])
    print("\tNumber of COMPLEX: ", table[table['variant_type'] == 'COMPLEX'].shape[0])
    print("\tNumber of INS SV: ", SVtable[SVtable['variant_type'] == 'INS'].shape[0])
    print("\tNumber of DEL SV: ", SVtable[SVtable['variant_type'] == 'DEL'].shape[0])
    print("\tNumber of COMPLEX SV: ", SVtable[SVtable['variant_type'] == 'COMPLEX'].shape[0])
    print("\tNumber of SNV: ", table[table['variant_type'] == 'SNV'].shape[0])

def filtered_variant_stats(table: pandas.DataFrame):
    bub_ids = set(table['bub_id'].tolist())
    print("\tNumber of variants filtered out: ", table.shape[0])
    print("\tNumber of bubbles involved: ", len(bub_ids))
    alt_lengths = table[['bub_id', 'alt_allele_len']]
    alt_lengths = alt_lengths.drop_duplicates()
    n_alleles = []
    for alt in alt_lengths['alt_allele_len'].tolist():
        n_alleles.append(len(alt.split(',')))
    
    print(pandas.Series(n_alleles).describe())
    pass

parser = argparse.ArgumentParser(prog='plot-vcf-stats.py', description="Making plots out of the VCF callset stats.")
parser.add_argument('-table', metavar='TABLE', help='statistics table')
parser.add_argument('-output', metavar='OUTPUT', help='output directory')
args = parser.parse_args()

table = pandas.read_csv(args.table, sep='\t', header=0)
print('Count before filtering:')
print_basic_stats(table)

print("\nFiltering variant records without any samples genotyped.")
filtered_variant_stats(table[table['all_total_alleles'] == 0])

table = table[table['all_total_alleles'] != 0]
print('\nCount after filtering:')
print_basic_stats(table)

# creating Table object for plotting
stats = Table(table)

# plotting HWE for all samples for all variant types
plot_all_hwe(stats._stats, args.output+'/hwe.png')
# plotting AF plot for all samples for all variant types
plot_all_panelAF_vs_callsetAF(stats._stats, args.output+'/af.png')

# plotting HWE and AF for pop-wise samples for all variant types
for pop in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
    plot_pop_hwe(stats._stats, args.output+'/hwe.%s.png'%(pop), pop=pop)
    plot_pop_panelAF_vs_callsetAF(stats._stats, args.output+'/af.%s.png'%(pop), pop=pop)

# making a dataframe with only SVs
stats_SV = stats.subset(is_SV=True)
# plotting HWE for all samples for SVs only
plot_all_hwe(stats_SV, args.output+'/hwe.SVonly.png')
# plotting AF plot for all samples for all variant types
plot_all_panelAF_vs_callsetAF(stats_SV, args.output+'/af.SVonly.png')

# plotting HWE and AF for pop-wise samples for SVs only
for pop in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
    plot_pop_hwe(stats_SV, args.output+'/hwe.%s.SVonly.png'%(pop), pop=pop)
    plot_pop_panelAF_vs_callsetAF(stats_SV, args.output+'/af.%s.SVonly.png'%(pop), pop=pop)


# making a dataframe with different vtypes
# SV only
stats_SV_vtype = {}
stats_SV_vtype['INS'] = stats.subset(is_SV=True, var_type='INS')
stats_SV_vtype['DEL'] = stats.subset(is_SV=True, var_type='DEL')
stats_SV_vtype['COMPLEX'] = stats.subset(is_SV=True, var_type='COMPLEX')
# all variants
stats_vtype = {}
stats_vtype['INS'] = stats.subset(var_type='INS')
stats_vtype['DEL'] = stats.subset(var_type='DEL')
stats_vtype['COMPLEX'] = stats.subset(var_type='COMPLEX')

# plotting HWE and AF
for vtype in ['INS', 'DEL', 'COMPLEX']:
    plot_all_hwe(stats_vtype[vtype], args.output+'/hwe.%s.png'%(vtype))
    plot_all_panelAF_vs_callsetAF(stats_vtype[vtype], args.output+'/af.%s.png'%(vtype))
    plot_all_hwe(stats_SV_vtype[vtype], args.output+'/hwe.%s.SVonly.png'%(vtype))
    plot_all_panelAF_vs_callsetAF(stats_SV_vtype[vtype], args.output+'/af.%s.SVonly.png'%(vtype))
