# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 09:35:35 2016

@author: dan
"""

from capacity_usage import CAPACITY_USAGE
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
def despine(ax, fontsize=15):
    ax.tick_params(right=0, top=0, direction='out', labelsize=fontsize)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel(ax.get_xlabel(), size=15)
    ax.set_ylabel(ax.get_ylabel(), size=15)
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    

flux = pd.DataFrame.from_csv("../data/mmol_gCDW_h.csv")
abundance =  pd.DataFrame.from_csv("../data/g_gCDW.csv")
cu = CAPACITY_USAGE(flux, abundance)

#%%
plt.figure(figsize=(6,6))
ax = plt.axes()
plt.scatter(cu.gr, cu.MCU)
despine(ax)
ax.set_xticks(np.arange(0,0.71,0.1))
ax.set_yticks(np.arange(0,0.71,0.1))
for c in cu.cs:
    ax.annotate(c, (cu.gr[c], cu.MCU[c]))
    
ax.set_xlabel("growth rate [$h^{-1}$]", size=15)
ax.set_ylabel("MCPU", size=15)
plt.tight_layout()
plt.savefig("../figures/png/figure5.png", dpi=200)