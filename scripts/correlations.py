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
ECU.dropna(thresh=3, inplace=True)
for c in conditions:
    y = ECU[c].copy().dropna()
    x = E_by_reac[c].copy().loc[y.index]    
    plt.scatter(x,y,c='cadetblue',edgecolor='')
    plt.xscale('log')

plt.xlabel('enzyme level [$mg \cdot gCDW^{-1}$]', size=15)
plt.ylabel('enzyme capacity usage', size=15)
[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]
#
#plt.axis('square')
#plt.xlim([1e-5,1e2])
#plt.ylim([-0.02,1])
#plt.annotate('R$^2$ = %.02f'%sr[0]**2, (0.2,0.85), xycoords='figure fraction',size=15)
#plt.tight_layout()
#plt.show()

#%%
#fig, ax = plt.subplots(1, 1)
#g = sb.jointplot(x, y, kind="kde", color="m", ax=ax)
#g.plot_joint(plt.scatter, c="w", s=30, linewidth=1, marker="+")
#g.ax_joint.collections[0].set_alpha(0)
#g.set_axis_labels("$X$", "$Y$");
#ax.set_xscale("log")

#sb.jointplot(np.log10(x), y, kind="kde")
#plt.close()
'''
#E to efficiency#
bg = '0.95'
#plt.figure(figsize=(8,8))
ax = plt.axes(axisbg=bg)
x = E_by_reac.dropna()
y = efficiency.dropna(thresh=3)
x = x.loc[y.index]

x['reactions']=x.index
x = pd.melt(x, id_vars='reactions')
x.columns = ['reactions','condition','E']

y['reactions']=y.index
y = pd.melt(y, id_vars='reactions')
y.columns = ['reactions','condition','CPU']

xy = pd.merge(x, y, on=['reactions', 'condition'])

gr.index.name = 'condition'
xy = pd.merge(xy, pd.DataFrame(gr), left_on = ['condition'],right_index=True)

#logx = np.log(x)
#logy = np.log(y)
#pr, sr = pearsonr(logx,y), spearmanr(logx,y)
xy.sort(columns=['growth rate [h-1]'], inplace=True)
xy.boxplot(column='CPU', by=['growth rate [h-1]'])

sb.jointplot(np.log10(xy["E"]), xy["CPU"], kind="kde")

#plt.xlabel('enzyme level [$mg \cdot gCDW^{-1}$]', size=15)
#plt.ylabel('enzyme capacity usage', size=15)
#[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
#[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]

#plt.axis('square')
#ax.set_xlim(0,1)
#plt.ylim([0,1])
#plt.annotate('R$^2$ = %.02f'%sr[0]**2, (0.2,0.85), xycoords='figure fraction',size=15)
#plt.tight_layout()
#plt.savefig("../res/Figure_S1.svg")
'''