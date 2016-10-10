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


fig = plt.figure(figsize=(6,6))
ax = plt.axes()

xy = df.loc[:,['reaction', 'Tm Protein','kmax', '# genes in reaction']].dropna(how='any')
xy = xy[xy['# genes in reaction']==1]
xy.drop_duplicates(subset='reaction', inplace=True)
xy.set_index('reaction', inplace=True)
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