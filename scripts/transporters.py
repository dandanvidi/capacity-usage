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

expression_CV = pd.DataFrame.from_csv(ppath+"supporting_data/ecoli_Schmidt_et_al_2015_CV.csv")[conds]
expression_CV.replace(np.nan,0, inplace=True)
expression_CV = expression_CV / 100 
expression_std = expression_CV.mul(mg_gCDW.mean(axis=1),axis=0)

umol_gCDW_min = get_umol_gCDW_min_from_pFVA(pFVA)
umol_gCDW_min = umol_gCDW_min.T[conds]

E = mg_gCDW
V = umol_gCDW_min
SA = specific_actitivy(V,E,model)
E_by_reac = V/SA

carbon = 'xyl_D'
c = gc[gc.media_key==carbon].index
transport = [r.id for r in model.reactions if carbon+'_p' in map(str,r.reactants)]
x0 = V.loc[transport][c]
y0 = E_by_reac.loc[transport][c]
#uptake_reactions = 