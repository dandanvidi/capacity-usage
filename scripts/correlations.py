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
from scipy.stats import pearsonr, spearmanr

gc = gc[gc.reference == 'Schmidt et al. 2015']
gc = gc[gc.strain == 'BW25113']
conditions = gc.dropna(subset=['media_key']).index & pFBA.columns
gr = gc['growth rate [h-1]'][conditions]
gr.sort()
conditions = gr.index

copies_fL = proteomics[conditions]
mmol_gCDW_h = pFBA[conditions]
#
mg_gCDW = convert_copies_fL_to_mg_gCDW(copies_fL)
umol_gCDW_min = mmol_gCDW_h * 1000 / 60

E = mg_gCDW
V = umol_gCDW_min


SA = specific_activity(V,E,model)

ECU = enzyme_capacity_usage(SA)
E_by_reac = (umol_gCDW_min/SA).loc[SA.index]

#E to efficiency#
bg = '0.95'
plt.figure(figsize=(8,8))
ax = plt.axes(axisbg=bg)
x = ECU.copy()
y = E_by_reac[x.columns].loc[x.index]
x = x.stack().values
y = y.stack().values

plt.scatter(y,x,c='cadetblue',edgecolor='')
plt.xscale('log')

plt.xlabel('enzyme level [$mg \cdot gCDW^{-1}$]', size=15)
plt.ylabel('enzyme capacity usage', size=15)
[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]
plt.xlim([1e-5,1e1])
plt.ylim([-0.02,1.02])
plt.tight_layout()
#plt.show()

#%%
for i,r in enumerate(E_by_reac.index):
    ry =  SA.loc[r]
    rx = V.loc[r]    
    R = pearsonr(rx,ry)
    if R[0]>0.3:
        plt.figure()
        plt.scatter(rx/rx.mean(),ry/ry.mean())
        plt.title(r)
        
    
#%%