# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 13:18:51 2016

@author: dan
"""
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.path.expanduser('~/git/kvivo_max/scripts/'))
from catalytic_rates import rates
import pandas as pd
import numpy as np
from helper import *

R = rates()
model = R.model
reactions = [r.id for r in model.reactions]
gc = R.gc[R.gc.media_key>0]
index = pd.MultiIndex.from_product([gc.index, ['maximum', 'minimum']])
fluxes = pd.DataFrame(index=index, columns=reactions)

for i,c in enumerate(gc.iterrows()):
    model = R.model
    gr = c[1]['growth rate [h-1]']
    cs = c[1]['media_key']
    ur = c[1]['uptake rate [mmol gCDW-1 h-1]']
    if np.isnan(ur):
        ur = 18.5
    fva = perform_pFVA(model, cs, gr, ur, reactions,fraction_of_optimum=1.00)
    fva = pd.DataFrame.from_dict(fva)
    fluxes.loc[c[0],'maximum'] = fva.loc['maximum']
    fluxes.loc[c[0],'minimum'] = fva.loc['minimum']
#    fluxes.index = pd.MultiIndex.from_product([c[0], fluxes.index])

fluxes.to_csv('../data/flux_variability_[mmol_gCDW_s]_relaxation=0.csv')

