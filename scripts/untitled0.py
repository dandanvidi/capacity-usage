# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 09:51:21 2015

@author: dan
"""
import pandas as pd
from helper import *
from scipy.optimize import curve_fit
#proteomics = metadata()
conditions = pd.DataFrame.from_csv('../data/growth_conditions.csv')
conditions = conditions[conditions['reference']=='Schmidt et al. 2015']
conditions = conditions[conditions['strain']=='BW25113']
copies_fL = pd.DataFrame.from_csv('../data/meta_abundance[copies_fL].csv')
copies_fL = copies_fL[conditions.index]

ribosome_genes = genes_by_function('Ribosome')
morethanone = {'b3985':2, 'b3986':4}
ribosome_stoichiometry = {g:morethanone[g] if g in morethanone else 1 for g in ribosome_genes}
ribosome_stoichiometry = pd.DataFrame.from_dict(ribosome_stoichiometry.items()).set_index(0)
del ribosome_stoichiometry.index.name
ribosome_stoichiometry = ribosome_stoichiometry[1]
ribosomal_genes_abundance = copies_fL.loc[ribosome_genes] # copies_fL
ribosomal_genes_weighted_abundance = ribosomal_genes_abundance.mul(ribosome_stoichiometry,axis=0).dropna()
ribosome_complex_abundance = (ribosomal_genes_weighted_abundance.sum() / 
                             ribosome_stoichiometry[ribosomal_genes_weighted_abundance.index].sum()) #copies_fL
                             

x = conditions['growth rate [h-1]']
y = ribosome_complex_abundance/1000

b_to_length = {row[0:5]:int(row[75:84]) for row in open("../data/all_ecoli_genes.txt", 'r')}
gene_length = pd.DataFrame.from_dict(b_to_length.items()).set_index(0)
gene_length = gene_length[1]
aa_abundance = copies_fL.mul(gene_length,axis=0).sum() #aa_fL
v = aa_abundance * x / 3600 #aa_fL_s
k_robosome = v/ribosome_complex_abundance

plt.figure(figsize=(6,6))
fontsize=15
ax = plt.axes()
plt.scatter(x, y,s=50,c='r',edgecolor='')
popt, pcov = curve_fit(lambda a,b,x: a*x+b, x, y)
intercept, slope = popt
plt.plot(x, x*slope+intercept)
ax.set_xlim(0)
ax.set_ylim(0)

ax.set_xlabel(r'growth rate $\left[ h^{-1} \right]$', size=fontsize)
ax.set_ylabel(r'ribosome abudance $\left[ \frac{copies}{fL} \times 10^3\right]$', size=fontsize)

[tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]
plt.tight_layout()
plt.savefig('../res/ribosome abundance by growth rate.pdf')

