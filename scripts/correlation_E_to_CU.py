# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 13:39:20 2016

@author: dan
"""

import pandas as pd
from capacity_usage import CAPACITY_USAGE
import matplotlib.pyplot as plt
#import seaborn as sns 
import sys, os
import numpy as np
from scipy.stats import ranksums, wilcoxon, pearsonr, spearmanr

cmap = plt.cm.Blues

flux = pd.DataFrame.from_csv("../data/mmol_gCDW_h.csv")
abundance =  pd.DataFrame.from_csv("../data/g_gCDW.csv")
cu = CAPACITY_USAGE(flux, abundance, shared_reactions=False)

x0 = cu.CU
y0 = cu.E

#%%

#'''plot mean values'''
#
#x1 = x0.mean(axis=1)
##x1.replace([np.nan, np.inf], 0, inplace=True)
#zeros1 = x1[x1==0]
#ones1 = x1[x1==1]
#x = x1[x1>0]
#x = x[x<1]
#y1 = y0.mean(axis=1)
#y = y1.loc[x.index]
#x.name = 'CU'
#y.name = 'E'
#xy = pd.concat([x,y], axis=1)
#
#xy.sort_values('CU',axis=0, inplace=True)
#
##sns.despine()
#plt.figure()
#ax = plt.axes()
#plt.scatter(x, y,edgecolor='', c='dodgerblue')
#
#plt.plot(xy['CU'], xy['E'].rolling(window=20, center=True).mean(), 'r')
#plt.ylim(1e-7, 1e2)
#plt.xlim(-.02, 1.02)
#plt.yscale('log')
#ax.spines['right'].set_color('none')
#ax.spines['top'].set_color('none')
#ax.spines['top'].set_color('none')
#ax.tick_params(direction='out')
#ax.xaxis.tick_bottom()
#ax.yaxis.tick_left()
#
#ax.scatter(zeros1, y1[zeros1.index], c='0.8', edgecolor='', zorder=0)
#ax.scatter(ones1, y1[ones1.index], c='0.8', edgecolor='', zorder=0)
#
#fs = 15
##plt.legend(loc=2, fontsize=fs)
#
#[tick.label.set_fontsize(fs) for tick in ax.xaxis.get_major_ticks()]
#[tick.label.set_fontsize(fs) for tick in ax.yaxis.get_major_ticks()]
#
#ax.set_xlabel('average capacity usage', fontsize=fs)
#ax.set_ylabel(r'average E $\left[\frac{mg}{gCDW}\right]$', fontsize=fs)
#ax.set_yticks([1e-6, 1e-4, 1e-2, 1e0, 1e2])
#plt.tight_layout()
#plt.savefig('E_to_CU_corr_mean.png')

#%%

'''plot all values'''

x2 = x0.stack(dropna=False)
#x2.replace([np.nan, np.inf], 0, inplace=True)
zeros2 = x2[x2==0]
ones2 = x2[x2==1]
x = x2
y2 = y0.stack(dropna=True)

y = y2.loc[x.index]
x.name = 'CU'
y.name = 'E'
xy = pd.concat([x,y], axis=1)

xy.sort_values('CU',axis=0, inplace=True)

#sns.despine()
plt.figure()
ax = plt.axes()
plt.scatter(x, y,edgecolor='', c='dodgerblue')

plt.plot(xy['CU'], xy['E'].rolling(window=350, center=True).median(), 'r')
plt.ylim(1e-7, 1e2)
plt.xlim(-.02, 1.02)
plt.yscale('log')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['top'].set_color('none')
ax.tick_params(direction='out')
ax.xaxis.tick_bottom()
ax.yaxis.tick_left()

ax.scatter(zeros2, y2[zeros2.index], c='0.8', edgecolor='', zorder=0)
ax.scatter(ones2, y2[ones2.index], c='0.8', edgecolor='', zorder=0)

fs = 15
#plt.legend(loc=2, fontsize=fs)

[tick.label.set_fontsize(fs) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fs) for tick in ax.yaxis.get_major_ticks()]

ax.set_xlabel('capacity usage', fontsize=fs)
ax.set_ylabel(r'E $\left[\frac{mg}{gCDW}\right]$', fontsize=fs)
ax.set_yticks([1e-6, 1e-4, 1e-2, 1e0, 1e2])
plt.tight_layout()
plt.savefig('../figures/png/E_to_CU_corr_all.png', dpi=200)
