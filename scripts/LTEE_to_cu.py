# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 14:09:58 2016

@author: dan
"""

import re, csv
import pandas as pd
from collections import defaultdict
from capacity_usage import CAPACITY_USAGE

flux = pd.DataFrame.from_csv("../data/mmol_gCDW_h.csv")
abundance =  pd.DataFrame.from_csv("../data/g_gCDW.csv")
cu = CAPACITY_USAGE(flux, abundance)

lenski = pd.DataFrame.from_csv("../data/Lenski_genes_[Pelosi_etal_2006].csv")
new_index = [re.split(r"[()]", s)[1] for s in lenski.index]
new_Q1 = [float(re.split(r" ", s)[0]) for s in lenski.Quantitatione1]
new_Q2 = [float(re.split(r" ", s)[0]) for s in lenski.Quantitatione2]
lenski.index = new_index
lenski.Quantitatione1 = new_Q1
lenski.Quantitatione2 = new_Q2
lenski = lenski[~lenski.index.duplicated(keep='first')]

#lenski.loc[x, 'Quantitatione1'] = 1 / lenski.loc[x, 'Quantitatione1']
#lenski.loc[x, 'Quantitatione2'] = 1 / lenski.loc[x, 'Quantitatione2']
#lenski.loc[:,'lower in'] = ['Evol']*len(lenski)
#%%
for g in lenski.index:
    with open('../data/all_ecoli_genes.txt', 'rb') as csvfile:
        f = csv.reader(csvfile)
        for i, s in enumerate(f):
            if g in s[0]:
                lenski.loc[g, 'bnumber'] = s[0][:5]
                break
#%%
                
out = pd.DataFrame(index=cu.rxns.keys(), columns=lenski.columns)
for g in lenski.index:
    b = lenski.bnumber[g]
    try:
        b = cu.model.genes.get_by_id(b)
        reactions = map(str, b.reactions)
        for r in reactions:
            out.loc[r] = lenski.loc[g]
    except KeyError:
        continue
out.dropna(subset=['Quantitatione1'], inplace=True)
out = out[['Quantitatione1', 'Quantitatione2', 'lower in']]
x = out.loc[:,['Quantitatione1', 'Quantitatione2']]
out.loc[:, 'Quantitatione'] = x.mean(axis=1)
out = out[['lower in', 'Quantitatione']]
out.loc[:, 'CU on glucose'] = (1-cu.CU.loc[out.index, 'glucose']) * cu.E['glucose'] + 1e-3
out.loc[:, 'Evol/Anc'] = out.Quantitatione
for x in out.index:
    if out.loc[x, 'lower in'] == 'Anc':
        if out.loc[x, 'Evol/Anc'] < 1:
            out.loc[x, 'Evol/Anc'] = 1 / out.loc[x, 'Evol/Anc']
    if out.loc[x, 'lower in'] == 'Evo':
        if out.loc[x, 'Evol/Anc'] > 1:
            out.loc[x, 'Evol/Anc'] = 1 / out.loc[x, 'Evol/Anc']
            
#%%
import matplotlib.pyplot as plt
plt.scatter(out['CU on glucose'], out['Evol/Anc'], edgecolor='')
plt.axhline(y=1, color='r')
#plt.xlim(-0.1,1.1)
plt.ylim(-0.1,6)
plt.xscale('log')
plt.ylabel('Evolved /Ancestor')
plt.xlabel('capacity usage')
out.replace(np.inf, np.nan, inplace=True)
x = out[out['CU on glucose']<0.1].index
print out.loc[x].mean()
cu.E.loc['MALS']
x = out[out['Evol/Anc']<1]
print len(x)