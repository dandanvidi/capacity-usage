# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 17:15:01 2016

@author: dan
"""

import re, pulp
import pandas as pd
import matplotlib.pyplot as plt
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
#from ..scripts.despine_axes import despine

def despine(ax, fontsize=15):
    ax.tick_params(right=0, top=0, direction='out', labelsize=fontsize)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel(ax.get_xlabel(), size=15)
    ax.set_ylabel(ax.get_ylabel(), size=15)
     
#%%
rid_mapping = pd.DataFrame.from_csv("../source/rid_mapping_cobra_2_Gerosa.csv")
MFA = pd.DataFrame.from_csv('../source/mmol_gCDW_hr_[Gerosa et al 2015].csv',
                                index_col=1)

MFA_std = pd.DataFrame.from_csv('../source/mmol_gCDW_hr_stdev_[Gerosa et al 2015].csv', 
                                    index_col=1)

conditions = pd.DataFrame.from_csv("../data/conditions.csv")
conditions = conditions[conditions.media_key>0]
conditions.sort_values('growth rate Gerosa [h-1]', inplace=True)
cs = conditions.index

#%%

measured_flux = pd.DataFrame(columns=cs, index=rid_mapping.index)
measured_flux_stdev = pd.DataFrame(columns=cs, index=rid_mapping.index)

for row in MFA.iterrows():
    if not re.findall("[+-]", row[0]):
        for r in row[0].split(';'):
            cobra_reactions = rid_mapping[rid_mapping['gerosa_reaction_id']==r]
            for r_cobra in cobra_reactions.index:
                v = row[1]
                measured_flux.loc[r_cobra] = v
                measured_flux_stdev.loc[r_cobra] = MFA_std.loc[row[0]]

measured_flux.dropna(inplace=True)
measured_flux_stdev.dropna(inplace=True)

#%%
model = create_cobra_model_from_sbml_file('../source/iJO1366.xml')
all_reactions = map(str, model.reactions)
all_metabolites = map(str, model.metabolites)

mmol_gCDW_h = pd.DataFrame(columns=cs, index=measured_flux.index)
for c in cs:
    cobra_c = conditions.loc[c, 'media_key']    
    gr = conditions.loc[c, 'growth rate Gerosa [h-1]']
    flux_meas = measured_flux[c]
    flux_stderr = measured_flux_stdev[c]

    # load fresh copy of model    
    model = create_cobra_model_from_sbml_file('../source/iJO1366.xml')
    
    # redefine sole carbon source uptake reaction in mmol/gr/h
    model.reactions.get_by_id('EX_glc_e').lower_bound = 0 
    model.reactions.get_by_id('EX_' + cobra_c + '_e').lower_bound = -1000 
            
    # set growth rate according to measurements
    biomass = "Ec_biomass_iJO1366_WT_53p95M"
    growth_rate = model.reactions.get_by_id(biomass)
    growth_rate.upper_bound = gr
    growth_rate.lower_bound = gr
    
    bounds_df = pd.DataFrame(index=all_reactions,columns=['lb','ub'])
    m = model.to_array_based_model()
    bounds_df.loc[all_reactions, 'lb'] = m.lower_bounds
    bounds_df.loc[all_reactions, 'ub'] = m.upper_bounds
    
    # initialize LP problem
    pulp_solver = pulp.CPLEX(msg=0)
    lp = pulp.LpProblem("MOMA", pulp.LpMinimize)
    
    v_pred = pulp.LpVariable.dicts('v_pred', all_reactions)
    v_meas = pulp.LpVariable.dicts('v_meas', all_reactions)
    v_resid = pulp.LpVariable.dicts('residual', all_reactions)   
    
    # add flux bounds
    for i in all_reactions:
        lp += (v_pred[i] >= bounds_df.loc[i, 'lb']), 'lower_bound_%s' % i
        lp += (v_pred[i] <= bounds_df.loc[i, 'ub']), 'upper_bound_%s' % i
        
    # add constraint for each measured reaction i:
    # |v_meas[i] - flux_meas[i]| <= flux_stderr[i]
    #  v_resid[i] >= |v_pred[i] - v_meas[i]|
    for i in flux_meas.index:
        lp += (v_meas[i] <= flux_meas[i] + flux_stderr[i]), 'measured_upper_%s' % i
        lp += (v_meas[i] >= flux_meas[i] - flux_stderr[i]), 'measured_lower_%s' % i
        lp += (v_pred[i] - v_resid[i] <= v_meas[i]), 'abs_diff_upper_%s' % i
        lp += (-v_pred[i] - v_resid[i] <= -v_meas[i]), 'abs_diff_lower_%s' % i

    # Some reactions in Gerosa et al. 2015 share constraints with other reactions
    # here we manually constrain their fluxes according to measuremnts. 

    # Acetate exchange
    lp += (v_meas['ACt2rpp'] + v_meas['ACS'] <= MFA.loc['PTAr+ACS', c] + MFA_std.loc['PTAr+ACS', c])
    lp += (v_meas['ACt2rpp'] + v_meas['ACS'] >= MFA.loc['PTAr+ACS', c] - MFA_std.loc['PTAr+ACS', c])
    
    # PFK/FBP reversible reaction
    lp += (v_meas['PFK'] - v_meas['FBP'] <= MFA.loc['PFK-FBP', c] + MFA_std.loc['PFK-FBP', c])
    lp += (v_meas['PFK'] - v_meas['FBP'] >= MFA.loc['PFK-FBP', c] - MFA_std.loc['PFK-FBP', c])
    
    # MDH/MQO alternative
    lp += (v_meas['MDH'] + v_meas['MDH2'] <= MFA.loc['MDH+MQO', c] + MFA_std.loc['MDH+MQO', c])
    lp += (v_meas['MDH'] + v_meas['MDH2'] >= MFA.loc['MDH+MQO', c] - MFA_std.loc['MDH+MQO', c])
    
    # ME alternative
    lp += (v_meas['ME1'] + v_meas['ME2'] <= MFA.loc['ME1+ME2', c] + MFA_std.loc['ME1+ME2', c])
    lp += (v_meas['ME1'] + v_meas['ME2'] >= MFA.loc['ME1+ME2', c] - MFA_std.loc['ME1+ME2', c])
      
    # set the objective to minimize sum_i abs_diff[i]
    objective = pulp.lpSum(v_resid.values())
    lp.setObjective(objective)        
        
    # add stoichiometric constraints for all internal metabolites: S_int * v = 0
    for i,j in enumerate(m.S):
        row = [l * v_pred[all_reactions[k]] for k,l in zip(j.rows[0],j.data[0])]
        lp += (pulp.lpSum(row) == 0), 'mass_balance_%s' % all_metabolites[i]        
    
    lp.solve()       

    # append fluxes to new dataframe
    MEAS_FLUX_L = 'measured fluxes from Gerosa et al.'
    MEAS_STDEV_L = 'standard deviation'
    PRED_FLUX_L = 'projected fluxes'
    RESID_L = 'residual'
    
    fluxes_df = pd.DataFrame(index=all_reactions)
    fluxes_df.loc[flux_meas.index, MEAS_FLUX_L] = flux_meas
    fluxes_df.loc[flux_meas.index, MEAS_STDEV_L] = flux_stderr
    
    fluxes_df.loc[all_reactions, PRED_FLUX_L] = \
                    map(lambda i: pulp.value(v_pred[i]), all_reactions)
    fluxes_df.loc[measured_flux.index, RESID_L] = \
                    map(lambda i: pulp.value(v_resid[i]), measured_flux.index)
    
    mmol_gCDW_h[c] = fluxes_df.loc[measured_flux.index, PRED_FLUX_L]

    #%%
    # normalize all fluxes to the biomass flux (i.e. set it to 1)
    fluxes_df /= pulp.value(v_pred[biomass])
    fig = plt.figure(figsize=(6,6))
    ax = plt.axes()
    
    fluxes_df.plot(kind='scatter', x=MEAS_FLUX_L, y=PRED_FLUX_L, 
                   xerr=MEAS_STDEV_L, ax=ax, linewidth=0, s=20,
                   color=(0.7,0.2,0.5))
    xlim, ylim = (ax.get_ylim(), ax.get_ylim())
    plt.axis('equal')
    plt.plot(xlim, ylim)
    plt.xlim(xlim)
    plt.ylim(ylim)
    despine(ax)

    ax.set_title(c, size=15)
    for i in flux_meas.index:
        xy = fluxes_df.loc[i, [MEAS_FLUX_L, PRED_FLUX_L]]
        if fluxes_df.loc[i, RESID_L] > 2:
            ax.annotate(i, xy,
                        fontsize=10, color='darkslategrey')
    fig.savefig('../res/flux_projections/flux_projection_on_%s.pdf' %c)

mmol_gCDW_h.to_csv('../data/flux projections[mmol_gCDW_h].csv')
