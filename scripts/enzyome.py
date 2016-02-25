import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.path.expanduser('~/git/kvivo_max/scripts/'))
#from catalytic_rates import rates
import pandas as pd
import numpy as np
from helper import *
import uncertainties.unumpy as unumpy 
from scipy.stats import pearsonr

gc = gc[gc.reference == 'Schmidt et al. 2015']
gc = gc[gc.strain == 'BW25113']
conditions = gc.dropna(subset=['media_key']).index & pFBA.columns
gr = gc['growth rate [h-1]'][conditions]
gr.sort()
conditions = gr.index

copies_fL = proteomics[conditions]
mmol_gCDW_h = pFBA[conditions]

mg_gCDW = convert_copies_fL_to_mg_gCDW(copies_fL)
umol_gCDW_min = mmol_gCDW_h * 1000 / 60

#use FVA instead of pFBA
#umol_gCDW_min = get_umol_gCDW_min_from_pFVA(pFVA)
#umol_gCDW_min = umol_gCDW_min.T[conditions]


SA = specific_activity(umol_gCDW_min,mg_gCDW,model)
E_by_reac = (umol_gCDW_min/SA).loc[SA.index]

ECU = enzyme_capacity_usage(SA)
MC = metabolic_capacity(umol_gCDW_min,mg_gCDW,model)
MCU = metabolic_capacity_usage(umol_gCDW_min,mg_gCDW,model)[conditions]

plt.scatter(gr,MCU)
r,p = pearsonr(gr,MCU)
plt.plot(range(0,2), range(0,2))
ax = plt.axes()
ax.set_xlim(0,0.8)
ax.set_ylim(0,0.8)
plt.annotate('$r^2 = %.2f$'%r**2,(0.15,0.875),xycoords='figure fraction',size=15)
ax.set_xlabel('growth rate [h$^{-1}$]', size=15)
ax.set_ylabel('capacity usage of metabolism', size=15)
[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]
#
#
#
#
#
#
#
#
#
#expression_CV = pd.DataFrame.from_csv(ppath+"supporting_data/ecoli_Schmidt_et_al_2015_CV.csv")[conds]
#expression_CV.replace(np.nan,0, inplace=True)
#expression_CV = expression_CV / 100 
#expression_std = expression_CV.mul(mg_gCDW.mean(axis=1),axis=0)
#
##umol_gCDW_min = get_umol_gCDW_min_from_pFVA(pFVA)
#umol_gCDW_min = pFBA[conds] * 1000 / 60
##umol_gCDW_min = umol_gCDW_min.T[conds]
#
#E = mg_gCDW
#V = umol_gCDW_min
#
#SA = specific_actitivy(V,E,model)
#E_by_reac = V/SA
#capacity = get_metabolic_capacity(V,E,model)
#usage = get_usage(V,E,model)
#capacity_usage = get_capacity_usage(V,E,model)
#
##standard_error = bootstrap_capacity_usage_error(V,E,model,iterations=10)
#
#bg = '0.95'
#fig = plt.figure(figsize=(8,8))
#ax = plt.axes(axisbg=bg)
#for i, c in enumerate(conds):
#    color='#ff4d4d'
##    if gc['growth mode'][c] == 'batch':
##        ax.annotate(gc['media_key'][c],(gr[c],capacity_usage[c]+0.01),
##                    ha='center',va='baseline',size=15)
##    elif gc['growth mode'][c] == 'chemostat':
##        color = '0.5'
#    plt.scatter(gr[c],capacity_usage[c],c=color,s=80,edgecolor='none')
##    ax.errorbar(gr[c],capacity_usage[c],standard_error[c],c='k')
#
#
#ax.set_xlim(0,0.8)
#ax.set_ylim(0,0.8)

#
#
#Vmax = E_by_reac.mul(SA.max(axis=1), axis=0)
#
#c = gc.index[0]
#x = []
#for c in gr.index:
#    GR_max = get_maximal_growth_rate(model, Vmax, c)
#    x.append(GR_max)
#    
#x = np.array(x)
##ax.hlines(capacity_usage[gr.index], gr, x, lw=8, color='y', zorder=0)
#
#plt.tight_layout()
#
#
#
##plt.savefig('../res/Figure_2.svg')
