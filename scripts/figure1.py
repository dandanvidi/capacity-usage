# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 16:00:14 2016

@author: dan
"""
import seaborn
from itertools import combinations, product
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.path.expanduser('~/git/kvivo_max/scripts/'))
#from catalytic_rates import rates
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import pandas as pd
import numpy as np
from helper import *
import uncertainties.unumpy as unumpy 
import matplotlib.cm as cm

model = create_cobra_model_from_sbml_file('../data/iJO1366.xml')
convert_to_irreversible(model)

rxns = {r.id:r for r in model.reactions}

gc = pd.DataFrame.from_csv('../data/growth_conditions.csv')
gc = gc[gc.media_key>0]
gc = gc[gc.reference == 'Schmidt et al. 2015']
gc = gc[gc.strain == 'BW25113']

pFVA = pd.DataFrame.from_csv("../data/flux_variability_[mmol_gCDW_h].csv", header=[0,1]).T
gr = gc['growth rate [h-1]'][pFVA.index.levels[0] & gc.index]
gr.sort()
conds = gr.index
ppath = "../../proteomics-collection/"

copies_fL = pd.DataFrame.from_csv(ppath+"meta_abundance[copies_fL].csv")[conds]
mg_gCDW = convert_copies_fL_to_mg_gCDW(copies_fL)

umol_gCDW_min = get_umol_gCDW_min_from_pFVA(pFVA)
umol_gCDW_min = umol_gCDW_min.T[conds]

V = umol_gCDW_min
SA = specific_actitivy(V,mg_gCDW,model)
E = V.div(SA).astype(float)
E = E.div(E.sum())

combs = [(c1, c2) for (c1, c2) in combinations(conds, 2)]
dlogEV = pd.DataFrame(index=SA.index, columns=combs)
median_dlogEV = pd.DataFrame(index=conds, columns=conds)
for x in combs:
    c1, c2 = x
    dlogEV[x] = np.log(E[c2] / E[c1])  /  np.log(V[c2] / V[c1])
    tmp = dlogEV[x].copy()
    tmp.replace([-np.inf, np.inf], np.nan, inplace=True)
    tmp.dropna(inplace=True)
    median_dlogEV[c2][c1] = tmp.median()
    
median_dlogEV = median_dlogEV.astype(float)
mask = np.zeros_like(median_dlogEV)
mask[np.triu_indices_from(mask)] = True

plt.figure()
ax=plt.axes()
labels = range(len(conds))
seaborn.heatmap(median_dlogEV, cmap='RdYlGn_r',vmin=-1, vmax=1,linewidths=2.5,ax=ax,
                xticklabels=labels,yticklabels=labels, square=True)
plt.yticks(rotation=0) 
plt.xticks(rotation=90) 
#plt.savefig("../res/figure1_red_green.svg")
plt.close()

def draw_foldchange_hist(ax, array, c):
    ax.hist(array, histtype='stepfilled', 
            bins=np.arange(-1, np.max(array) + 0.1, 0.25),                                                         
            color='#CCCCB2', edgecolor='none')
    ax.axvline(array.median(), 0, 1, c=c, zorder=10, ls='-', lw=4)
    

x2 = ('ACE_BATCH_mu=0.3_S','GLC_BATCH_mu=0.58_S')
x1 = ('GLC_CHEM_mu=0.12_S','GLC_CHEM_mu=0.50_S')
hists = [x1,x2]
colors = [[1.,0.9686274509803922,0.6980392156862745,1.],[0.91372549,0.33333333,0.21960784,1.]]

f, axes = plt.subplots(2, 1, sharey=True)
for i, x in enumerate(hists):    
    tmp = dlogEV[x].copy()
    tmp = tmp[tmp>-5]
    tmp.replace([-np.inf, np.inf], np.nan, inplace=True)
    tmp.dropna(inplace=True)    
    draw_foldchange_hist(axes[i], tmp, colors[i])
#    seaborn.distplot(tmp,ax=ax, bins=40)
    axes[i].set_xlim(-1,2)
    axes[i].set_ylim(0,160)
    print tmp.median()
    
f.set_figheight(4)
f.set_figwidth(3)
plt.tight_layout() 
#plt.savefig("../res/histograms.svg")

#matric = [dE_to_dV[x].dropna().median() for x in combs]