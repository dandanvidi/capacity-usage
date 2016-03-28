import matplotlib.pyplot as plt
import sys, os
import pandas as pd
import numpy as np
from helper import *
import uncertainties.unumpy as unumpy 
from scipy.stats import pearsonr
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import seaborn
from collections import Counter

MFA = pd.DataFrame.from_csv('../data/MFA_Gerosa_et_al_2015.csv', header=0, index_col=0)

conditions = MFA.columns & mg_gCDW.columns

#gc = gc[gc.reference == 'Schmidt et al. 2015']
#gc = gc[gc.strain == 'BW25113']
#gc = gc[gc['growth mode'] == 'batch']
#x = gc.copy().dropna(subset=['comments']).index
#gc.drop(x, inplace=True)
#conditions = gc.dropna(subset=['media_key']).index & pFBA.columns
gr = gc['growth rate [h-1]'][conditions]
gr.sort()
conditions = gr.index

copies_fL = proteomics[conditions]
mmol_gCDW_h = pFBA[conditions]
#
mg_gCDW = convert_copies_fL_to_mg_gCDW(copies_fL)
umol_gCDW_min = mmol_gCDW_h * 1000 / 60
MCU = metabolic_capacity_usage(umol_gCDW_min,mg_gCDW,model)[conditions]

#MCU.sort()
plt.figure()
ax = plt.axes()
plt.bar(range(len(conditions)), MCU*100)
ax.set_xticks(np.arange(len(conditions))+0.4)
ax.set_xticklabels([gc["media_key"][c] for c in MCU.index], rotation='vertical')
#plt.ylim(30, 100)
plt.xlim(0,len(conditions))
plt.ylabel("capacity usage [%]", size=15)
[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]

plt.tight_layout()