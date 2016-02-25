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

#expression_CV = pd.DataFrame.from_csv(ppath+"supporting_data/ecoli_Schmidt_et_al_2015_CV.csv")[conds]
#expression_CV.replace(np.nan,0, inplace=True)
#expression_CV = expression_CV / 100 
#expression_std = expression_CV.mul(mg_gCDW.mean(axis=1),axis=0)

umol_gCDW_min = get_umol_gCDW_min_from_pFVA(pFVA)
umol_gCDW_min = umol_gCDW_min.T[conds]

E = mg_gCDW
V = umol_gCDW_min


SA = specific_actitivy(V,E,model)
efficiency = get_efficiency(V,E,model)
E_by_reac = V/SA

#E to efficiency#
bg = '0.95'
plt.figure(figsize=(8,8))
ax = plt.axes(axisbg=bg)
x = E_by_reac.median(axis=1).dropna()
y = efficiency.dropna(thresh=3)
y = y.median(axis=1)[x.index].dropna()
x = x.loc[y.index]
logx = np.log(x)
logy = np.log(y)
pr, sr = pearsonr(logx,y), spearmanr(logx,y)
plt.scatter(x,y,c='cadetblue',edgecolor='')
plt.xscale('log')

plt.xlabel('enzyme level [$mg \cdot gCDW^{-1}$]', size=15)
plt.ylabel('enzyme capacity usage', size=15)
[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]

plt.axis('square')
plt.xlim([1e-5,1e2])
plt.ylim([-0.02,1])
plt.annotate('R$^2$ = %.02f'%sr[0]**2, (0.2,0.85), xycoords='figure fraction',size=15)
plt.tight_layout()
plt.show()

#%%
fig, ax = plt.subplots(1, 1)
g = sb.jointplot(x, y, kind="kde", color="m", ax=ax)
g.plot_joint(plt.scatter, c="w", s=30, linewidth=1, marker="+")
g.ax_joint.collections[0].set_alpha(0)
g.set_axis_labels("$X$", "$Y$");
ax.set_xscale("log")

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