# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 13:29:23 2016

implementation of FBA with molecular crowding - FBAwMC

Beg, Q. K. et al. Intracellular crowding defines the mode and sequence of 
substrate uptake by Escherichia coli and constrains its metabolic activity. 
Proc. Natl. Acad. Sci. U. S. A. 104, 12663–12668 (2007).

@author: dan
"""

from cobra.core import Metabolite, Reaction
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import pandas as pd
import numpy as np
from helper import *

#ppath = "../../proteomics-collection/"
pFVA = pd.DataFrame.from_csv("../data/flux_variability_[mmol_gCDW_h].csv", header=[0,1]).T
pFBA = pd.DataFrame.from_csv("../data/flux[mmol_gCDW_h].csv")
gc = pd.DataFrame.from_csv('../data/growth_conditions.csv')
gc = gc[gc.media_key>0]
gc = gc[gc.reference == 'Schmidt et al. 2015']
gc = gc[gc.strain == 'BW25113']
gr = gc['growth rate [h-1]'][pFVA.index.levels[0] & gc.index]
gr.sort()
conds = gr.index

copies_fL = pd.DataFrame.from_csv("../data/meta_abundance[copies_fL].csv")[conds]
mg_gCDW = convert_copies_fL_to_mg_gCDW(copies_fL)
#umol_gCDW_min = get_umol_gCDW_min_from_pFVA(pFVA)
#umol_gCDW_min = umol_gCDW_min.T[conds]

umol_gCDW_min = pFBA[conds] * 1000 / 60


E = mg_gCDW
V = umol_gCDW_min

model = create_cobra_model_from_sbml_file('../data/iJO1366.xml')
convert_to_irreversible(model)
rxns = {r.id:r for r in model.reactions}
rxns['Ec_biomass_iJO1366_core_53p95M'].objective_coefficient = 0            
rxns['Ec_biomass_iJO1366_WT_53p95M'].objective_coefficient = 1      
rxns['EX_glc_e'].lower_bound = -1000
crowding = Metabolite(id='crowding')
model.add_metabolites(crowding)                

M = pd.Series(get_complex_molecular_weight(model)) / 1000 # kDa
SA = specific_actitivy(V,E,model)
kmax = SA.max(axis=1) # umol/mg/min
kmax = convert_SA_to_kcat(kmax, M) #1/s

ECU = get_efficiency(V,E,model)
ECU = ECU.loc[rxns.keys()][conds]
ECU[ECU==0] = np.nan
ECU_pad = ECU.replace({c:{np.nan:ECU[c].mean()} for c in conds}).copy()
ECU_rand = get_rand_ECU(ECU_pad,model)

kmax_pad = kmax.copy()
kmax_pad[kmax_pad==0] = np.nan
kmax_pad = kmax_pad.replace(np.nan, kmax.median()) # assign median kmax for reactions without kmax
#rand_ECU[rand_ECU>0] = 1.
#rand_ECU[rand_ECU>=0] = .6

v_specific = 0.73 # mL/g
density = 0.34 # g/mL
v_molar = M * v_specific * 1000 # mL/mol
a = 1 / ECU_pad.mul(kmax_pad,axis=0)
a = density * a.mul(v_molar,axis=0)
a = a / (3600*1000)

c = conds[-4]
for r in model.reactions:
    r.add_metabolites({crowding:a[c][r.id]}) 

crowding_counter = Reaction(name='crowding_counter')
crowding_counter.add_metabolites(metabolites={crowding:-1})                
model.add_reaction(crowding_counter) 
crowding_counter.upper_bound = 1 * 0.5 * 0.44 *6.05 / 21.7# molecular crowding constraint

model.optimize()
print (0.7 / model.solution.f) * 60

ACU = ECU.mean()
plt.figure(figsize=(6,6))
ax = plt.axes()
plt.scatter(gr, ACU)
x = gr.copy()
x[0] = 0
x[1] = 1
plt.plot(x,x,'r')
plt.xlim(0,0.9)
plt.ylim(0,0.9)
plt.show()
plt.close()

#
x = pd.Series(index=map(str,model.reactions), data=model.solution.x)
y = pFBA[c]
#y[y<1e-10] = np.nan
#y.dropna(inplace=True)
x = x[y.index]
plt.figure(figsize=(6,6))
ax = plt.axes()
z = (y/x).dropna()
plt.scatter(range(len(z)), z)
#z = y/x
z = z.replace(np.inf, np.nan).dropna()
z.sort()


#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlim(1e-6, 1e4)
#ax.set_ylim(1e-6, 1e4)
plt.tight_layout()
#for r in model.reactions:
''' 



CU = np.arange(0.01,1, 0.1)
gr = []
aiz = []
#for cu in CU:
cu = 0.6
model = create_cobra_model_from_sbml_file('../data/iJO1366.xml')
convert_to_irreversible(model)
rxns = {r.id:r for r in model.reactions}

rxns['Ec_biomass_iJO1366_core_53p95M'].objective_coefficient = 0            
rxns['Ec_biomass_iJO1366_WT_53p95M'].objective_coefficient = 1      
rxns['EX_glc_e'].lower_bound = -1000
for r in model.reactions:
    M_i = E_weights[r.id] # g/mol
    try:
        kmax_i = convert_SA_to_kcat(kmax[r.id], M_i/1000)
    except KeyError:
        kmax_i = 10
    
    v_i_molar = M_i * v_i_specific # mL/mol
    a_i = density * v_i_molar / (kmax_i * cu)
    a_i  = a_i / 3600 / 1000
    r.add_metabolites({crowding:a_i})    

#        model.remove_reactions([r])



model.optimize()
gr.append(model.solution.f)
aiz.append(model.solution.x_dict['crowding_counter'])
gr = model.solution.f
print (0.7/gr) * 60

#plt.figure()
#plt.plot(CU,gr,'r')
#plt.xlabel('capacity usage')
#plt.ylabel('growth rate')
#plt.figure()
#plt.plot(aiz,CU,'b')
'''