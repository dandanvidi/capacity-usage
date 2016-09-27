# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 13:43:14 2016

@author: dan
"""

import pandas as pd
from capacity_usage import CAPACITY_USAGE
import matplotlib.pyplot as plt
#import seaborn as sns 
import sys, os
import numpy as np
from scipy.stats import ranksums, wilcoxon
from gauge import gauge, degree_range

def despine(ax, fontsize=15):
    ax.tick_params(right=0, top=0, direction='out', labelsize=fontsize)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel(ax.get_xlabel(), size=15)
    ax.set_ylabel(ax.get_ylabel(), size=15)
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    
flux = pd.DataFrame.from_csv('../data/Gerosa et al [mmol_gCDW_h].csv')
flux_sd = pd.DataFrame.from_csv('../data/Gerosa et al stdev[mmol_gCDW_h].csv')
abundance =  pd.DataFrame.from_csv("../data/g_gCDW.csv")
abundance_cv = pd.DataFrame.from_csv('../data/g_CDW_CV.csv')
cu = CAPACITY_USAGE(flux, abundance, shared_reactions=False)

#%%
r = 'G6PDH2r'
b = list(cu.model.reactions.get_by_id(r).genes)[0].id
cp = cu.CU.loc[r]
sa = cu.SA.loc[r]
E  = cu.E.loc[r]
Emin = cu.Emin.loc[r]
E_sd = E*abundance_cv.loc[b]/100
E_sd = E_sd[E.index]
v  = cu.umol_gCDW_min.loc[r]
vmax  = cu.vmax.loc[r]
v_sd = flux_sd.loc[r]*1000/60
mu = cu.gr

#%%

fig, (ax1, ax2) = plt.subplots(2,1, figsize=(12,4), sharex=True)
despine(ax1)
despine(ax2)

ax1.bar(np.arange(len(v)), v.values, yerr=v_sd, width=0.3, color='#786721', 
        edgecolor='', ecolor='k')
ax1.bar(np.arange(len(vmax))-0.3/1.5, vmax.values, width=0.3, color='#E9DDAF',
        zorder=0)

ax2.bar(np.arange(len(E))-0.3/1.5, E.values, yerr=E_sd, width=0.3, color='#577061',
        edgecolor='', ecolor='k')
ax2.bar(np.arange(len(Emin)), Emin.values, width=0.3, color='0.8', edgecolor='')
#E.plot(kind='bar', ax=ax2, yerr=E_sd, width=0.3, color='#577061')
#Emin.plot(kind='bar', ax=ax2, width=0.3, color='0.8', zorder=10)

#ax1.set_ylim(0,7)
#ax1.set_yticks(np.arange(0,8,3.5))
ax1.set_ylim(0,200)
ax2.set_ylim(0,.200)


ax2.set_yticks(np.arange(0,.201,.0500))
ax1.set_ylabel(r'v [$\frac{mmol}{gCDW \cdot hr}$]')
ax2.set_ylabel(r'E [$\frac{mg}{gCDW}$]')

ax2.set_xticklabels(['']+list(cu.cs)+[''])
plt.tight_layout(h_pad=3)

plt.savefig('../res/zwf_bars.svg')
#ax.set_ylim(0,0.05)


#plt.figure()
#ax = plt.axes()

#ax.set_ylim(0,3)



#%%
#
#N=30
#conds = cu.cs
#for p in conds:
#    gauge(N=N, colors='viridis', 
#          cat=cp[p], top_title=p, 
#          title=r'CU=%.2f'%cp[p],
#          fname='../res/gauge/zwf_%s.svg'%p)          