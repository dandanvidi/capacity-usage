# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 10:27:12 2016

@author: dan
"""

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
cw = gc['cell width [um]'][pFVA.index.levels[0] & gc.index]
cl = gc['cell length [um]'][pFVA.index.levels[0] & gc.index]
s2v = get_surface_to_volume_ratio(cl,cw)[2].dropna()
gr = gr[s2v.index]
gr.sort()

plt.figure()
ax = plt.axes()
plt.scatter(gr,s2v,c='coral',edgecolor='',s=50)
ax.set_xlim(0,1)
ax.set_ylim(0,10)
ax.set_ylabel('surface to volume [$\mu$m$^{-1}$]', size=15)
ax.set_xlabel('growth rate [h$^{-1}$]', size=15)
[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]
plt.savefig("../res/s2v is constant across growth rates.svg")