# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 13:20:57 2016

@author: dan
"""

import pandas as pd
import numpy as np
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from cobra.manipulation.modify import convert_to_irreversible

conditions = pd.DataFrame.from_csv("../data/conditions.csv")
conditions = conditions[conditions.media_key>0]
conditions.sort_values('growth rate Gerosa [h-1]', inplace=True)
cs = conditions.index

projections = pd.DataFrame.from_csv('../data/flux projections[mmol_gCDW_h].csv')

model = create_cobra_model_from_sbml_file('../source/iJO1366.xml')
convert_to_irreversible(model)
mmol_gCDW = pd.DataFrame(index = map(str, model.reactions),
                         columns = cs)
for c in projections.columns:

    cobra_c = conditions.loc[c, 'media_key']    
    gr = conditions.loc[c, 'growth rate Gerosa [h-1]']
    model = create_cobra_model_from_sbml_file('../source/iJO1366.xml')
    convert_to_irreversible(model)

    # redefine sole carbon source uptake reaction in mmol/gr/h
    model.reactions.get_by_id('EX_glc_e').lower_bound = 0 
    model.reactions.get_by_id('EX_' + cobra_c + '_e').lower_bound = -1000 
            
    # set growth rate according to measurements
    biomass = "Ec_biomass_iJO1366_WT_53p95M"
    model.change_objective(biomass)
    growth_rate = model.reactions.get_by_id(biomass)
    growth_rate.upper_bound = gr
#    growth_rate.lower_bound = gr
    
    revs = []   
    for rid in projections.index:
        v = projections.loc[rid, c]
        r_f = model.reactions.get_by_id(rid)
        if  v < 0:
            r_b = model.reactions.get_by_id(rid+'_reverse')
            revs.append(r_b)
            r_f.upper_bound = 0
            r_f.lower_bound = 0
            r_b.upper_bound = -v
            r_b.lower_bound = -v
        else:
            r_f.upper_bound = v
            r_f.lower_bound = v

    f = optimize_minimal_flux(model, already_irreversible=True)
    mmol_gCDW[c].update(pd.Series(f.x_dict))

umol_gCDW_min = mmol_gCDW * 1000 / 60
mmol_gCDW.to_csv('../data/mmol_gCDW_hr.csv')
umol_gCDW_min.to_csv('../data/umol_gCDW_min.csv')