# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 10:48:23 2016

@author: dan
"""
from sklearn.cluster import KMeans
from capacity_usage import CAPACITY_USAGE
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import ranksums

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
cu = CAPACITY_USAGE(flux, abundance, shared_reactions=False)

biosyn_mean = []
ccm_mean = []
ranksum_pvalues = {}
dists = []
for c in cu.cs:
    x = cu.CU[c]
    x.replace(np.inf, np.nan, inplace=True)
    x.dropna(inplace=True)
    subsystems = pd.Series(index=x.index, data=[cu.rxns[r].subsystem for r in x.index])
    for k,v in cu.get_master_groups().iteritems():
        subsystems.replace({i:k for i in v}, inplace=True)
        
    s = set(subsystems.values)
    l = []
    for j, s0 in enumerate(['central metabolism', 'biosynthesis']):
        i = subsystems[subsystems==s0].index
        l.append(x[i])
        if s0 == 'central metabolism':
            ccm_mean.append(x[i].median())
        if s0 == 'biosynthesis':
            biosyn_mean.append(x[i].median())
        dists.append(x[i])
    ranksum_pvalues[c] = ranksums(l[0], l[1])[1]


#%%
#ind = np.arange(len(cu.gr))  # the x locations for the groups
#width = 0.35       # the width of the bars
#
#fig, ax = plt.subplots()
#rects1 = ax.bar(ind, ccm_mean, width, color='r')
#
#rects2 = ax.bar(ind + width, biosyn_mean, width, color='y')
#
## add some text for labels, title and axes ticks
##ax.set_ylabel('Scores')
##ax.set_title('Scores by group and gender')
#ax.set_xticks(ind + width)
#ax.set_xticklabels(cu.cs, rotation=90, fontsize=15)
#despine(ax)
#ax.legend((rects1[0], rects2[0]), ('Men', 'Women'))
plt.figure(figsize=(8,6))
ax = plt.axes()
#for i, l0 in enumerate(dists): 
#    print i
x = [dists[m] for m in range(len(dists))]
positions=[[j,j+0.3] for j in range(len(dists)/2)]
tmp = []
for i in positions:
    tmp += i

#boxprops = dict(linestyle='--', linewidth=3, color='darkgoldenrod')
box = plt.boxplot([dists[i] for i in range(len(dists))], 
            positions=tmp,
            patch_artist=True,
#            notch=True,
#            boxprops=boxprops,
            widths=0.3)

colors = ['lightblue', 'tan']*7
labels= ['Central carbon', 'Biosynthesis']
for i, (patch, color) in enumerate(zip(box['boxes'], colors)):
    patch.set_facecolor(color)
    patch.set_edgecolor('k')
ax.set_xticks([np.mean(positions[n]) for n in range(len(positions))])
ax.set_xticklabels(cu.cs, size=15, rotation=90)
despine(ax)
ax.set_ylim(0,1)
plt.savefig("../res/biosynthesis_vs_ccm.svg")
#x = [(dists[i][0], dists[i][1]) for i in range(len(dists))]
#ax.boxplot([dists[1][0], dists[1][1]], widths=0.05, positions=[1,1+0.1])

