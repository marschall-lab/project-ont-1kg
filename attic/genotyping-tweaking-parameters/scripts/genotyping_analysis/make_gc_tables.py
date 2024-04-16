import sys
import subprocess
import pandas as pd
from collections import defaultdict
import pickle
import matplotlib.pyplot as plt

path = '/home/samarendra/mount/hpc-medbioinf/users/spani/results/1000GP/genotype-tweaking/genotype_concordance/'
ls_command = 'ls '+path
process = subprocess.Popen(ls_command.split(), stdout=subprocess.PIPE)
ls_out, ls_err = process.communicate()

files = [f for f in ls_out.decode('utf-8').rstrip().split('\n') if ".pkl" in f]

wgc = defaultdict(lambda: defaultdict(lambda: None))
typed = defaultdict(lambda: defaultdict(lambda: None))
prec = defaultdict(lambda: defaultdict(lambda: None))
recall = defaultdict(lambda: defaultdict(lambda: None))
for file in files:
    f = file.split(".")[0].split("_")
    ws = f[1][2:]
    wh_flag = f[2]
    vtype = f[3]
    data = None
    with open(path+file, 'rb') as tmp:
        data = pickle.load(tmp)
    wgc[vtype][(ws, wh_flag)] = data['weighted_concordance']
    typed[vtype][(ws, wh_flag)] = data['typed']
    prec[vtype][(ws, wh_flag)] = data['precision']
    recall[vtype][(ws, wh_flag)] = data['recall']


for vtype, values in wgc.items():
    wgc_table = {'snp': [], 'indel': []}
    typed_table = {'snp': [], 'indel': []}
    prec_table = {'snp': [], 'indel': []}
    recall_table = {'snp': [], 'indel': []}
    ws = [int(i[0]) for i in list(values.keys())]
    ws = sorted(list(set(ws)))
    for f in ['snp','indel']:
        for w in ws:
            wgc_table[f].append(round(wgc[vtype][(str(w),f)],5))
            typed_table[f].append(round(typed[vtype][(str(w),f)],5))
            prec_table[f].append(round(prec[vtype][(str(w),f)],5))
            recall_table[f].append(round(recall[vtype][(str(w),f)],5))
    if vtype in ['snp', 'all']:
        columns = ['snp', 'indel']
    else:
        columns = ['indel']
    rows = ws
    
    df = pd.DataFrame(wgc_table, index=pd.Index(rows), columns=pd.Index(columns))
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (4,2))
    ax = plt.subplot(111, frame_on=False) # no visible frame
    ax.xaxis.set_visible(False) # hide the x axis
    ax.yaxis.set_visible(False) # hide the y axis
    pd.plotting.table(ax, df, loc='center') # where df is your data frame
    plt.tight_layout()
    plt.savefig('/home/samarendra/mount/hpc-medbioinf/users/spani/results/1000GP/genotype-tweaking/genotype_concordance/tables/%s_wgc.png'%(vtype))
    plt.close()
    df = pd.DataFrame(typed_table, index=pd.Index(rows), columns=pd.Index(columns))
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (4,2))
    ax = plt.subplot(111, frame_on=False) # no visible frame
    ax.xaxis.set_visible(False) # hide the x axis
    ax.yaxis.set_visible(False) # hide the y axis
    pd.plotting.table(ax, df, loc='center') # where df is your data frame
    plt.tight_layout()
    plt.savefig('/home/samarendra/mount/hpc-medbioinf/users/spani/results/1000GP/genotype-tweaking/genotype_concordance/tables/%s_typed.png'%(vtype))
    plt.close()
    df = pd.DataFrame(prec_table, index=pd.Index(rows), columns=pd.Index(columns))
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (4,2))
    ax = plt.subplot(111, frame_on=False) # no visible frame
    ax.xaxis.set_visible(False) # hide the x axis
    ax.yaxis.set_visible(False) # hide the y axis
    pd.plotting.table(ax, df, loc='center') # where df is your data frame
    plt.tight_layout()
    plt.savefig('/home/samarendra/mount/hpc-medbioinf/users/spani/results/1000GP/genotype-tweaking/genotype_concordance/tables/%s_precision.png'%(vtype))
    plt.close()
    df = pd.DataFrame(recall_table, index=pd.Index(rows), columns=pd.Index(columns))
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (4,2))
    ax = plt.subplot(111, frame_on=False) # no visible frame
    ax.xaxis.set_visible(False) # hide the x axis
    ax.yaxis.set_visible(False) # hide the y axis
    pd.plotting.table(ax, df, loc='center') # where df is your data frame
    plt.tight_layout()
    plt.savefig('/home/samarendra/mount/hpc-medbioinf/users/spani/results/1000GP/genotype-tweaking/genotype_concordance/tables/%s_recall.png'%(vtype))
    plt.close()