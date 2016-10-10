# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 11:43:31 2016

@author: dan
"""

import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import pandas as pd
import seaborn

df = pd.DataFrame.from_csv("../cache/thermal_stability_metadata.csv")
    
fig = plt.figure(figsize=(6,6))
ax = plt.axes()

xy = df.loc[:,['bnumber','GLC_BATCH_mu=0.58_S', 'Protein Abundance']].dropna(how='any')
xy.drop_duplicates(subset='bnumber', inplace=True)
xy.set_index('bnumber', inplace=True)
x = xy.iloc[:,0]
y = xy.iloc[:,1]
ax.scatter(x, y, edgecolor='', c='#8080ff')
ax.set_xlim(1e-1, 1e5)
ax.set_ylim(1e-1, 1e5)
ax.set_xscale('log')
ax.set_yscale('log')

ax.annotate(r'$r^2=%.2f$'%spearmanr(x,y)[0]**2, (0.1,0.9), 
            xycoords='axes fraction', size=15)

ax.set_xlabel('Schmidt et al. [copies/fL]', size=15)
ax.set_ylabel('Leuenberger et al. [?]', size=15)
ax.tick_params(direction='out', labelsize=15)
plt.tight_layout()
plt.savefig("../res/ts_abundance.pdf")
#%%

fig = plt.figure(figsize=(6,6))
ax = plt.axes()

xy = df.loc[:,['reaction', 'Tm Protein','kapp_glucose', '# genes in reaction']].dropna(how='any')
xy = xy[xy['# genes in reaction']==1]
xy.drop_duplicates(subset='reaction', inplace=True)
xy.set_index('reaction', inplace=True)
x = xy.iloc[:,0]
y = xy.iloc[:,1]
ax.scatter(x, y, edgecolor='', c='#e6ac00')
ax.set_xlim(40, 70)
ax.set_ylim(1e-5, 1e5)
#ax.set_xscale('log')
ax.set_yscale('log')

ax.annotate(r'$r^2=%.2f$'%spearmanr(x,y)[0]**2, (0.1,0.9), 
            xycoords='axes fraction', size=15)

ax.set_xlabel('$T_m$ Protein'  , size=15)
ax.set_ylabel(r'$k_\mathrm{app}$', size=18)
ax.tick_params(direction='out', labelsize=15)
plt.tight_layout()
plt.savefig("../res/ts_tmvskapp.pdf")
#%%

fig = plt.figure(figsize=(6,6))
ax = plt.axes()

df.loc[:,'kcat/kapp'] = df['kcat'] / df['kapp_glucose']
xy = df.loc[:,['reaction', 'Tm Protein','kcat/kapp', '# genes in reaction']].dropna(how='any')
xy = xy[xy['# genes in reaction']==1]
xy.drop_duplicates(subset='reaction', inplace=True)
xy.set_index('reaction', inplace=True)
x = xy.iloc[:,0]
y = xy.iloc[:,1]
ax.scatter(x, y, edgecolor='')
ax.set_xlim(40, 70)
ax.set_ylim(1e-3, 1e3)
#ax.set_xscale('log')
ax.set_yscale('log')

ax.annotate(r'$r^2=%.2f$'%spearmanr(x,y)[0]**2, (0.1,0.9), 
            xycoords='axes fraction', size=15)

ax.set_xlabel('$T_m$ Protein', size=15)
ax.set_ylabel(r'$k_\mathrm{cat}$ / $k_\mathrm{app}$', size=18)
ax.tick_params(direction='out', labelsize=15)
plt.tight_layout()
plt.savefig("../res/ts_tmvskcat.pdf")

#%%
fig = plt.figure(figsize=(6,6))
ax = plt.axes()

df.loc[:,'kcat/kapp'] = df['kcat'] / df['kapp_glucose']
xy = df.loc[:,['reaction', 'Aggregation','kcat/kapp', '# genes in reaction']].dropna(how='any')
xy = xy[xy['# genes in reaction']==1]
xy.drop_duplicates(subset='reaction', inplace=True)
xy.set_index('reaction', inplace=True)
x = xy.iloc[:,0]
y = xy.iloc[:,1]
ax.scatter(x, y, edgecolor='', c='darkred')
ax.set_xlim(-1e-1, 1e-1)
ax.set_ylim(1e-3, 1e3)
#ax.set_xscale('log')
ax.set_yscale('log')

ax.annotate(r'$r^2=%.2f$'%spearmanr(x,y)[0]**2, (0.1,0.9), 
            xycoords='axes fraction', size=15)

ax.set_xlabel('Aggregation', size=15)
ax.set_ylabel(r'$k_\mathrm{cat}$ / $k_\mathrm{app}$', size=18)
ax.tick_params(direction='out', labelsize=15)
plt.tight_layout()
plt.savefig("../res/ts_aggregationvsresidual.pdf")
