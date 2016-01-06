# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 13:59:45 2015

@author: dan
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

fontsize = 15

gc = pd.DataFrame.from_csv("../data/growth_conditions.csv")
gc = gc[gc.reference=='Schmidt et al. 2015']
gr = gc['growth rate [h-1]'][gc.index]

data = pd.DataFrame.from_csv('../data/protein_abundance[mg_gCDW].csv', sep='\t')
conditions = data.columns & gr.index

mg_gCDW = data[conditions]



plt.figure(figsize=(8,5))
ax = plt.axes()
xstd = gc['growth rate std [h-1]'][conditions]
new = gc[gc.comments=='new'].index

for c in conditions:
    if c not in new:
        plt.scatter(gr[c], mg_gCDW.sum()[c],s=40,c='#d966ff',edgecolor='')
#        plt.errorbar(gr[c], mg_gCDW.sum()[c], xerr=xstd[c],c='k')
    else:
        print c
        plt.scatter(gr[c], mg_gCDW.sum()[c],s=40,edgecolor='',c='#ffcc33')
#        plt.errorbar(gr[c], mg_gCDW.sum()[c], xerr=xstd[c],c='m')
        
ax.set_xlabel(r'growth rate $\left[ h^{-1} \right]$', size=fontsize)
ax.set_ylabel(r'total protein $\left[ \frac{mg}{gCDW} \right]$', size=fontsize)



[tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]
    
ax.set_xlim(-0.1,2)
#ax.set_ylim(450,500)
#ax.grid()
plt.tight_layout()
plt.savefig('../res/proteome concentration by gr.png')

new_data = data[new]

plt.figure(figsize=(8,8))
ax = plt.axes()
plt.plot(data['GLC_BATCH_mu=0.58_S'],data['GLC_BATCH_mu=0.6_NEW'], 'ro', alpha=0.5)
print spearmanr(data['GLC_BATCH_mu=0.58_S'],data['GLC_BATCH_mu=0.6_NEW'])
logx = np.log(data['GLC_BATCH_mu=0.58_S'].replace(0,np.nan)).dropna()
logy = np.log(data['GLC_BATCH_mu=0.6_NEW'].replace(0,np.nan)).dropna()
shared = logx.index & logy.index
print pearsonr(logx.loc[shared],logy.loc[shared])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-7, 1e2)
ax.set_ylim(1e-7, 1e2)
[tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]
ax.set_xlabel(r'total protein $\left[ \frac{mg}{gCDW} \right]$', size=fontsize)
ax.set_ylabel(r'total protein $\left[ \frac{mg}{gCDW} \right]$ - new', size=fontsize)
plt.tight_layout()
plt.savefig('../res/comparing glucose.png')

plt.figure(figsize=(8,8))
ax = plt.axes()
plt.plot(data['FUM_BATCH_mu=0.42_S'],data['FUM_BATCH_mu=0.55_NEW'], 'co', alpha=0.5)

print spearmanr(data['FUM_BATCH_mu=0.42_S'],data['FUM_BATCH_mu=0.55_NEW'])

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-7, 1e2)
ax.set_ylim(1e-7, 1e2)
[tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]
ax.set_xlabel(r'total protein $\left[ \frac{mg}{gCDW} \right]$', size=fontsize)
ax.set_ylabel(r'total protein $\left[ \frac{mg}{gCDW} \right]$ - new', size=fontsize)
plt.savefig('../res/comparing fumarate.png')

