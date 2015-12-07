# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 13:59:45 2015

@author: dan
"""

import pandas as pd
import matplotlib.pyplot as plt

gc = pd.DataFrame.from_csv("../data/growth_conditions.csv")
gc = gc[gc.reference=='Schmidt et al. 2015']
gr = gc['growth rate [h-1]'][gc.index]

data = pd.DataFrame.from_csv('../data/protein_abundance[mg_gCDW].csv', sep='\t')
conditions = data.columns & gr.index

mg_gCDW = data[conditions]

plt.figure(figsize=(8,5))
ax = plt.axes()
plt.scatter(gr, mg_gCDW.sum(),s=40,edgecolor='',c='#01735C')
ax.set_xlabel(r'growth rate $\left[ h^{-1} \right]$', size=fontsize)
ax.set_ylabel(r'total protein $\left[ \frac{mg}{gCDW} \right]$', size=fontsize)

fontsize = 15

[tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]
    
ax.set_xlim(0,2)
ax.set_ylim(450,500)
ax.grid()

plt.savefig('../res/proteome concentration by gr.png')