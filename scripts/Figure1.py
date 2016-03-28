# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 15:13:32 2016

@author: dan
"""
from load_proteomics import PROTEOMICS
from capacity_usage import CAPACITY_USAGE
import pandas as pd

p = PROTEOMICS()
cu = CAPACITY_USAGE()

data = cu.ECU.dropna(how="any")
data.index.name = "reaction"

subsystems = pd.Series(index=data.index, data=[cu.rxns[r].subsystem for r in data.index])
subsystems.index.name = "reaction"
subsystems.name = "subsystem"


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
             
             
data = data.join(subsystems)

data.replace({"subsystem":{i:"amino acids biosynthesis" for i in AA_biosyn}}, inplace=True)
x = data.groupby("subsystem").mean()
#pivoted = data.pivot(index="reaction", columns="function", values=p.conditions.index)
#
#metabolism = p.gene_info[p.gene_info["Annotated functional COG class"] == "METABOLISM"]
#COG_metabolism = metabolism["Annotated functional COG group (description)"]
#metabolic_mass = p.mg_gCDW.loc[metabolism.index].sum() / p.mg_gCDW.sum() *100

'''
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
umol_gCDW_min[umol_gCDW_min<1e-10] = 0
#
SA = specific_activity(umol_gCDW_min,mg_gCDW,model)
E_by_reac = (umol_gCDW_min/SA).loc[SA.index]
y = E_by_reac.dropna(thresh=len(conditions)-2).copy()
x = umol_gCDW_min.loc[y.index].copy()


#for i in x.index:
i = "GLUDy_reverse"
i = 'DHQS'
logx, logy = np.log10(umol_gCDW_min.loc[i]), np.log10(E_by_reac.loc[i])
plt.figure()
plt.plot(logx, logy, 'ro')
plt.title(i)
#if pearsonr(logx, logy)[0]**2 > 0.75:        
#    plt.figure()
#    plt.plot(x.loc[i], y.loc[i], 'ro')
#    plt.title(i)
    
'''