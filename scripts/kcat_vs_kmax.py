# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 17:18:19 2016

@author: dan
"""

import pandas as pd
from capacity_usage import CAPACITY_USAGE
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.stats import pearsonr, spearmanr


flux = pd.DataFrame.from_csv("../data/mmol_gCDW_h.csv")
abundance =  pd.DataFrame.from_csv("../data/g_gCDW.csv")
cu = CAPACITY_USAGE(flux, abundance)
def configure_plot(ax, x_label='', y_label='', fontsize=15):
    ax.tick_params(right=0, top=0, direction='out', labelsize=fontsize)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel(ax.get_xlabel(), size=15)
    ax.set_ylabel(ax.get_ylabel(), size=15)
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    ax.set_xlabel(x_label, size=fontsize*1.3)
    ax.set_ylabel(y_label, size=fontsize*1.3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    

def add_labels(x, y, labels, ax, fig, fontsize=10, hide_overlap=True):
    ann = []
    for name in labels.index:
        if x[name]>y[name]:
            ann.append(ax.text(x[name], y[name]/1.1, labels[name], 
                                ha='center', va='top', zorder=5, size=fontsize))
        if x[name]<y[name]:
            ann.append(ax.text(x[name], y[name]*1.1, labels[name],
                                ha='center', va='bottom', zorder=5, size=fontsize))
                                    
        mask = np.zeros(fig.canvas.get_width_height(), bool)
        fig.canvas.draw()
        for i, a in enumerate(ann):
            bbox = a.get_window_extent()
            x0 = int(bbox.x0)
            x1 = int(math.ceil(bbox.x1))
            y0 = int(bbox.y0)
            y1 = int(math.ceil(bbox.y1))
        
            s = np.s_[x0:x1, y0:y1]
            if hide_overlap:
                if np.any(mask[s]):
                    a.set_visible(False)
                else:
                    mask[s] = True
            else:
                mask[s] = True 
#%%  
outliers = ['ORPT_reverse', 'PPKr_reverse', 'GLUTRR', 
                      'PGCD', 'FCLT', 'CPPPGO']

too_low_E = ['FCLT', 'GLUTRR', 'ORPT_reverse']
fig = plt.figure(figsize=(7,7))
ax = plt.axes()
lowerlim=1e-4
upperlim=1e5
y = cu.kmax_per_sec()
y.replace(0, np.nan, inplace=True)
reactions = cu.kcat.dropna().index & y.dropna().index - too_low_E #- outliers
x = cu.kcat[reactions]
y = y[reactions]
ax.scatter(x, y, color='#be9b7b')
ax.plot([lowerlim, upperlim], [lowerlim, upperlim], 'b')
configure_plot(ax, 
               x_label=r"$k_\mathrm{cat}\left[ s^{-1} \right]$",
               y_label=r"$k_\mathrm{max}^{vivo}\left[ s^{-1} \right]$")

labels = pd.DataFrame.from_csv('../data/kcat_data.csv').loc[reactions,'gene name']

add_labels(x,y,labels,ax, fig, hide_overlap=True)
ax.scatter(x[outliers], y[outliers], color='darkred')

plt.xlim(lowerlim,upperlim)
plt.ylim(lowerlim,upperlim)
plt.tight_layout()
plt.savefig("../res/kcat_vs_kmax_with_outliers.png", dpi=300)
print spearmanr(x,y)[0]**2, spearmanr(np.log(x),np.log(y))[0]**2

#%%
#out = pd.DataFrame(index=x.index)
#out.loc[:, 'kcat [umol/mg/min]'] = x
#out.loc[:, 'kmax [umol/mg/min]'] = y
#out.loc[:, 'kcat/kmax'] = x/y
#out.loc[:, 'log10[kcat/kmax]'] = np.log10(x/y)
#out.loc[:, 'function'] = [cu.rxns[r].name for r in out.index]
#out.index.name = 'reaction'
#out.to_csv("../res/kcat vs kmax.csv", sep='\t')