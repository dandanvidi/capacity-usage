import pandas as pd
import sys, os, csv
from scipy.optimize import curve_fit
from helper import *
from cobra.manipulation import remove_genes, undelete_model_genes
from cobra.flux_analysis import single_gene_deletion
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
sys.path.append(os.path.expanduser('~/git/kvivo_max/scripts/'))
from catalytic_rates import rates, perform_pFBA
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize


R = rates()
avogadro = 6.022*1e23
model = R.model
data = pd.DataFrame.from_csv('../data/protein_abundance[mg_gCDW].csv', sep='\t')
conditions = data.columns & R.gc.index

gr = R.gc['growth rate [h-1]'][conditions]

mg_gCDW = data[conditions]
protein_MW = data['MW[Da]']

mmol_gCDW = mg_gCDW.div(protein_MW,0)

ribosome_genes = genes_by_function('Ribosome')
morethanone = {'b3985':2, 'b3986':4}
ribosome_stoichiometry = {g:morethanone[g] if g in morethanone else 1 for g in ribosome_genes}
ribosome_stoichiometry = pd.DataFrame.from_dict(ribosome_stoichiometry.items()).set_index(0)
ribosome_stoichiometry.index.name = 'bnumber'
ribosome_stoichiometry = ribosome_stoichiometry[1]
ribosome_mmol_gCDW = mmol_gCDW.loc[ribo_genes]

ribosome_abundance = ribosome_mmol_gCDW.mul(ribosome_stoichiometry,axis=0).dropna()
ribosome_abundance = ribosome_abundance.sum() / ribosome_stoichiometry[ribosome_abundance.index].sum() #mmol_gCDW
ribosome_flux = (mg_gCDW / 110 * (gr + 1/20.) / 3600).sum()#mmol_AA_gCDW_s



ribo_kapp = (ribosome_flux / ribosome_abundance) / 0.95 #1/s


plt.figure()
a = gc[gc.comments=='new']
x = gr[gc.index - a.index]
y = ribosome_abundance[gc.index - a.index]
popt, pcov = curve_fit(lambda a,b,x: a*x+b, x, y)
intercept, slope = popt
plt.plot(gr, gr*slope+intercept)
plt.plot(gr, ribosome_abundance, 'o')

#res = minimize(lambda a, b: ribosome_flux - (a*ribosome_abundance - a*b), (intercept, 20))

#upper_limit_degredation_rate = (ribosome_abundance * 20 - ribosome_flux).max()#mmol_AA_gCDW_s
#upper_limit_degredation_rate = intercept * 20

ribo_kapp = ((ribosome_flux) / ribosome_abundance) / 0.95 #1/s




#ribo_mg_gCDW = mg_gCDW.loc[ribo_genes].dropna()
#rRNA_mg_gCDW = ribo_mg_gCDW.sum() * 2
##RNTPs_mmol_gCDW = rRNA_mg_gCDW / 340 * gr / 3600

#
#rnap_genes = {b for b,ko in b_to_KEGG.iteritems() if ko in rnap_KEGG}
#rnap_stoich = {'b3295':2}
#rnap_genes = rnap_stoich.keys()
#rnap_stoich = pd.DataFrame.from_dict(rnap_stoich.items()).set_index(0)
#rnap_stoich.index.name = 'bnumber'
#rnap_stoich = rnap_stoich[1]
#rnap_mg_gCDW = mg_gCDW.loc[rnap_genes].dropna()
#rnap_mmol_gCDW = mmol_gCDW.loc[rnap_genes].dropna()
#rnap_kapp = ribo_RNA_mmol_gCDW / (rnap_mmol_gCDW.mul(rnap_stoich,axis=0)).sum()/4
#
#
bremer = np.array([[0.6, 1, 1.5, 2, 2.5],[12,16, 18, 20, 21]])
gCDW_cell = np.array([148, 258, 433, 641, 865]) * 1e-15 #gCDW_cell
ribosomes_per_cell = 1e3 * np.array([6.8, 13.5, 26.3, 45.1, 72]) / avogadro * 1e3 # mmol_cell
mmol_ribosomes_per_gCDW = ribosomes_per_cell / gCDW_cell
bremer_kapp = (bremer[1] * mmol_ribosomes_per_gCDW + upper_limit_degredation_rate) / mmol_ribosomes_per_gCDW

plt.figure(figsize=(8,5))
fontsize = 15
ax = plt.axes()
ax.set_xlabel(r'growth rate $\left[ h^{-1} \right]$', size=fontsize)
ax.set_ylabel(r'ribosome rate $\left[ \frac{AA}{s} \right]$', size=fontsize)
#1
#
xstd = R.gc['growth rate std [h-1]'][conditions]
#

for c in conditions:
    try:
        np.isnan(R.gc['comments'][c])
        if R.gc['growth mode'][c] == 'batch':
            plt.scatter(gr[c], ribo_kapp[c],s=40,edgecolor='#668cff',c='',zorder=0)
            plt.errorbar(gr[c], ribo_kapp[c], xerr=xstd[c],c='#668cff')
        if R.gc['growth mode'][c] == 'chemostat':
            plt.scatter(gr[c], ribo_kapp[c],s=40,edgecolor='k',c='',zorder=1)
            plt.errorbar(gr[c], ribo_kapp[c], xerr=xstd[c],c='k')
    except:
        if R.gc['comments'][c] == 'new':
            plt.scatter(gr[c], ribo_kapp[c],s=40,edgecolor='r',c='',zorder=2)
            plt.errorbar(gr[c], ribo_kapp[c], xerr=xstd[c],c='r')
        else:
            plt.scatter(gr[c], ribo_kapp[c],s=40,edgecolor='r',c='',zorder=2)
            plt.errorbar(gr[c], ribo_kapp[c], xerr=xstd[c],c='r')

[tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]
    
#plt.scatter(bremer[0],bremer_kapp,edgecolor='y',c='',marker='D',s=40)
ax.set_xlim(.25,.7)
ax.set_ylim(-1,25)
plt.tight_layout()
plt.savefig('../res/ribosome rate by gr.pdf')


#unique = set([g.id for g in R.model.genes if len(g.reactions)==1])
#KO_set = {g for g in R.genes if g not in mg_gCDW.index}
#KO_set.remove('s0001')
#flux = R.flux_data[conditions]
#kapp = R.kapp[conditions]
#kmax = R.kmax['kmax per chain [s^-1]']
#saturation = kapp.div(kmax,axis=0)
#index = set(proteins.index) & set(map(lambda x:x.id,enzymes))
#index = set(mg_gCDW.index) & unique
#
#mass = pd.DataFrame(index=efficiency.index, columns=gc.index)
#
#for reac in mass.index:
#    r = R.rxns[reac]
#    genes = map(lambda x: x.id, r.genes)
#    try:
#        mass.loc[reac] = mg_gCDW.loc[genes].sum()
#    except:
#        continue
#mass.dropna(how='all', inplace=True)


#for c in conditions:
    
    
#    cond = R.gc.loc[c]
#    gr = R.gc['growth rate [h-1]'][c]
#    cs = R.gc['media_key'][c]
#    ur = R.gc['uptake rate [mmol gCDW-1 h-1]'][c]
#    flux = perform_pFBA(model, cs, gr, ur)
