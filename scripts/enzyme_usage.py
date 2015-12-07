# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 13:01:03 2015

@author: dan
"""

from catalytic_rates import rates
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

gc = pd.DataFrame.from_csv("../data/growth_conditions.csv")
gc = gc[gc.reference=='Schmidt et al. 2015']
gr = gc['growth rate [h-1]'][gc.index]

data = pd.DataFrame.from_csv('../data/protein_abundance[mg_gCDW].csv', sep='\t')
conditions = data.columns & gr.index
mg_gCDW = data[conditions]

R = rates()
kapp = R.kapp
kmax = R.kmax['kmax per chain [s^-1]']
#kcat = R.kcat['kcat per chain [s^-1]']

kmax_usage = kapp.div(kmax, axis=0).dropna(how='all')
kmax_usage = kmax_usage[gc.index & kmax_usage.columns]


effective_capacity = pd.DataFrame(index=kmax_usage.index, columns=conditions)
for reac in effective_capacity.index:
    r = R.rxns[reac]
    genes = map(lambda x: x.id, r.genes)
    try:
        effective_capacity.loc[reac] = mg_gCDW.loc[genes].sum() * kmax_usage.loc[reac]
    except:
        continue
        
effective_capacity.dropna(how='all', inplace=True)

glc_ec = effective_capacity['GLC_BATCH_mu=0.58_S'].dropna().astype(float)
ace_ec = effective_capacity['ACE_BATCH_mu=0.3_S'].dropna().astype(float)

from plot_types import cdf

plt.figure()
ax = plt.axes()
cdf(glc_ec, ax=ax)
cdf(ace_ec, ax=ax)
ax.set_xscale('log')
ax.set_xlim(1e-3,1e0)
#kcat_usage = kapp.div(kmax, axis=0).dropna(how='any')
#kcat_usage = kmax_usage[gc.index & kmax_usage.columns]