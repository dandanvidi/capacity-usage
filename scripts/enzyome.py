import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.path.expanduser('~/git/kvivo_max/scripts/'))
#from catalytic_rates import rates
import pandas as pd
import numpy as np
from helper import *
import uncertainties.unumpy as unumpy 
from scipy.stats import pearsonr
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import seaborn
from collections import Counter

conditions = gc.index
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
#
SA = specific_activity(umol_gCDW_min,mg_gCDW,model)
E_by_reac = (umol_gCDW_min/SA).loc[SA.index]
#
ECU = enzyme_capacity_usage(SA)[conditions]
SUBSYSTEMS = pd.Series(index=ECU.index, data=[rxns[rx].subsystem for rx in ECU.index])

#%%
from sklearn.cluster import KMeans

X = ECU.dropna(how='any').values

est = KMeans(n_clusters=11)
est.fit(X)
labels = est.labels_

reacs = list(ECU.dropna(how='any').index)
reacs.sort(key=dict(zip(reacs, labels)).get)

clustered = ECU.loc[reacs].copy()
conds = [gc['media_key'][c] for c in conditions]
bg = '0.95'
plt.figure(figsize=(10,10))
ax=plt.axes()
cmap = plt.cm.viridis
cmap.set_bad(bg,1.)
seaborn.heatmap(clustered, cmap=cmap,ax=ax,xticklabels=conds,yticklabels=False)

r2group = pd.DataFrame(index=reacs, columns=['group','class','value'])
r2group['class'] = labels
r2group['group'] = SUBSYSTEMS.loc[reacs].values
r2group['value'] = np.ones(len(reacs))
groups2class = pd.pivot_table(r2group, index='group', values='value', aggfunc='sum', columns='class').replace(np.nan, 0)

plt.show()
plt.close()
#for s in set(SUBSYSTEMS):
    
    
from scipy import stats
#oddsratio, pvalue = stats.fisher_exact(groups2class.values)

#plt.xticks(rotation=90)

#%%
MC = metabolic_capacity(umol_gCDW_min,mg_gCDW,model)
MCU = metabolic_capacity_usage(umol_gCDW_min,mg_gCDW,model)[conditions]

plt.figure(figsize=(7,7))
plt.scatter(gr,MCU)
r,p = pearsonr(gr,MCU)

#a = np.cov(gr, MCU)[0,0] / np.var(gr)
a, cov = curve_fit(lambda a,x:a*x, gr,MCU)
plt.plot([0, 1], [0, a*1])
ax = plt.axes()
ax.set_xlim(0,0.8)
ax.set_ylim(0,0.8)
plt.annotate('$r^2 = \, %.2f$\n$\mu_{max} = \, %.2f \, h^{-1}$'% (r**2, 1/a),
             (0.15,0.8),xycoords='figure fraction',size=15)
ax.set_xlabel('growth rate [h$^{-1}$]', size=15)
ax.set_ylabel('capacity usage of metabolism', size=15)
[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]

AA_biosyn = ['Alanine and Aspartate Metabolism',
             'Arginine and Proline Metabolism',
             'Cysteine Metabolism',
             'Glutamate Metabolism',
             'Glycine and Serine Metabolism',
             'Histidine Metabolism',
             'Methionine Metabolism',
             'Threonine and Lysine Metabolism',
             'Tyrosine, Tryptophan, and Phenylalanine Metabolism',
             'Valine, Leucine, and Isoleucine Metabolism']

central_carbon = ['Citric Acid Cycle',
                  'Glycolysis/Gluconeogenesis',
                  'Anaplerotic Reactions',
                  'Pentose Phosphate Pathway']
#%%
import colorsys

def rmsd(x, f):
    return ((f-x)**2).sum()
    
def RgbToHex(rgb):
	r, g, b = rgb
	r *= 255
	g *= 255
	b *= 255
	return '#%02x%02x%02x' % (r,g,b)
 
def ColorMap(items, saturation=0.7, value=0.95, hues=None):
	if hues is None:
		n = len(items)
		hues = np.arange(float(n)) / float(n)
	f = lambda h: RgbToHex(colorsys.hsv_to_rgb(h, saturation, value)) 
	rgbs = map(f, hues)
	return dict(zip(items, rgbs))
 


AAs = [r for r in SUBSYSTEMS.index if SUBSYSTEMS[r] in AA_biosyn]
CAR = [r for r in SUBSYSTEMS.index if SUBSYSTEMS[r] in central_carbon]
plt.plot(gr,E_by_reac.loc[AAs].sum())

#%%

BIOSYN_REACS = pd.DataFrame(index=E_by_reac.index,columns=['R2', 'slope'])
for r in BIOSYN_REACS.index:
    try:
        y = E_by_reac.loc[r][conditions]
        a, cov = curve_fit(lambda a,x:a*x, gr,y)
        BIOSYN_REACS['R2'][r] = pearsonr(gr,y)[0]**2
        BIOSYN_REACS['slope'][r] = a[0]
    except ValueError:
        continue
BIOSYN_REACS.dropna(how='any', inplace=True)
plt.hist(BIOSYN_REACS['R2'], bins=20)
#plt.figure()
#y = ECU.loc['CYSS'][conditions]
#plt.scatter(gr,y)
#a, cov = curve_fit(lambda a,x:a*x, gr,y)
#plt.plot([0,1],[0,1*a])
#l = f(a,gr)
#print pearsonr(gr,y)[0]**2, a


#%%
CU = ECU*E_by_reac
CU_sub = (CU[conditions].groupby(by=SUBSYSTEMS).sum() / 
         E_by_reac[conditions].groupby(by=SUBSYSTEMS).sum())
for s in set(SUBSYSTEMS):
    plt.figure(figsize=(10,10))
    plt.title(s, size=15)
    plt.ylabel('capacity usage [%]', size=15)
    plt.xlabel('growth rate [h-1]', size=15)
    plt.plot(gr, MCU, 'k*')
    plt.scatter(gr, CU_sub.loc[s])
    plt.xlim(0,1)
    plt.ylim(0,1)
    
    
#E_by_reac['subsystem'] = [rxns[rx].subsystem for rx in E_by_reac.index]
#E_by_reac = E_by_reac[E_by_reac.subsystem!='']
#ECU = ECU[ECU.subsystem!='']
#ECU_by_subsys = ECU.groupby('subsystem').sum()
#E_by_reac_by_subsys = E_by_reac.groupby('subsystem').sum()
#CU_by_subsys = CU_by_subsys[conditions]
#plt.figure()
#CU_by_subsys = ECU_by_subsys.div(E_by_reac_by_subsys)
#for s in CU_by_subsys.index:
#    plt.plot(gr,CU_by_subsys.loc[s],'r', marker='o')

#%%


#E_by_reac['subsystem'] = [rxns[rx].subsystem for rx in E_by_reac.index]
#E_by_reac = E_by_reac[E_by_reac.subsystem!='']
#ECU = ECU[ECU.subsystem!='']
#ECU_by_subsys = ECU.groupby('subsystem').sum()
#E_by_reac_by_subsys = E_by_reac.groupby('subsystem').sum()
#CU_by_subsys = CU_by_subsys[conditions]
#plt.figure()
#CU_by_subsys = ECU_by_subsys.div(E_by_reac_by_subsys)
#for s in CU_by_subsys.index:
#    plt.plot(gr,CU_by_subsys.loc[s],'r', marker='o')

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
