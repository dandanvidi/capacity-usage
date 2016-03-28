# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 11:19:09 2016

@author: dan
"""

import pandas as pd
import pulp
import os
from cobra.io.sbml import create_cobra_model_from_sbml_file
import matplotlib
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
import matplotlib.pyplot as plt
from cobra.manipulation.modify import convert_to_irreversible

matplotlib.rcParams['text.usetex'] = False

MEAS_FLUX_L = 'measured fluxes from Gerosa et al.'
MEAS_STDEV_L = 'standard deviation'
PRED_FLUX_L = 'projected fluxes'
RESID_L = 'residual'
gc = pd.DataFrame.from_csv("../data/carbon_sources.csv")
MFA = pd.DataFrame.from_csv('../data/MFA_Gerosa_et_al_2015.csv', header=0, index_col=0)
MFA_std = pd.DataFrame.from_csv('../data/MFA_Gerosa_et_al_2015_deviations.csv', header=0, index_col=0)


    #%%
if __name__ == '__main__':

    
    mmol_gCDW_h = {}    
    for condition in gc.index:
        mmol_gCDW_h[condition] = {}
        cs = gc["media_key"][condition] 
        gr = gc["growth rate [h-1]"][condition]
        gr_stdev = gc["growth rate stdev [h-1]"][condition]
        model = create_cobra_model_from_sbml_file('../data/iJO1366.xml')
        model.reactions.get_by_id('EX_glc_e').lower_bound = 0 
        model.reactions.get_by_id('EX_' + cs + '_e').lower_bound = -1000 # redefine sole carbon source uptake reaction in mmol/gr/h
        biomass = model.reactions.get_by_id("Ec_biomass_iJO1366_core_53p95M")
        biomass.upper_bound = gr
        biomass.lower_bound = gr
            
        all_reactions = map(str, model.reactions)
        all_metabolites = map(str, model.metabolites)
        m = model.to_array_based_model()
        bounds_df = pd.DataFrame(index=all_reactions,columns=['lower_bound','upper_bound'])
        bounds_df['lower_bound'] = m.lower_bounds
        bounds_df['upper_bound'] = m.upper_bounds

        #%%
        mfa = MFA[condition]        
        mfa_std = MFA_std[condition]
        
        df = pd.DataFrame(index=mfa.index, columns=['mean','stdev'])
        df['mean'] = mfa
        df['stdev'] = mfa_std

#%%
        flux_means = df.iloc[:, 0]
        flux_stderr = df.iloc[:, 1]
    
        fluxes_df = pd.DataFrame(index=all_reactions)
        fluxes_df.loc[flux_means.index, MEAS_FLUX_L] = flux_means
        fluxes_df.loc[flux_means.index, MEAS_STDEV_L] = flux_stderr

#%%
        pulp_solver = pulp.CPLEX(msg=0)
        lp = pulp.LpProblem("FLUX_L1", pulp.LpMinimize)
        
        measured_reactions = list(flux_means.index)
    
        v_pred = pulp.LpVariable.dicts('v_pred', all_reactions)
        v_meas = pulp.LpVariable.dicts('v_meas', measured_reactions)
        v_resid = pulp.LpVariable.dicts('residual', measured_reactions)   
        
        # add flux bounds
        for i in all_reactions:
            lp += (v_pred[i] >= bounds_df.loc[i, 'lower_bound']), 'lower_bound_%s' % i
            lp += (v_pred[i] <= bounds_df.loc[i, 'upper_bound']), 'upper_bound_%s' % i
            
        # add constraint for each measured reaction i:
        # |v_meas[i] - flux_means[i]| <= flux_stderr[i]
        # v_resid[i] >= |v_pred[i] - v_meas[i]|
        for i in measured_reactions:
            lp += (v_meas[i] <= flux_means[i] + flux_stderr[i]), 'measured_upper_%s' % i
            lp += (v_meas[i] >= flux_means[i] - flux_stderr[i]), 'measured_lower_%s' % i
            lp += (v_pred[i] - v_resid[i] <= v_meas[i]), 'abs_diff_upper_%s' % i
            lp += (-v_pred[i] - v_resid[i] <= -v_meas[i]), 'abs_diff_lower_%s' % i
            
        # also set the objective to be minimizing sum_i abs_diff[i]
        objective = pulp.lpSum(v_resid.values())
        lp.setObjective(objective)        
            
        # add stoichiometric constraints for all internal metabolites: S_int * v = 0
        for i,j in enumerate(m.S):
            row = [l * v_pred[all_reactions[k]] for k,l in zip(j.rows[0],j.data[0])]
            lp += (pulp.lpSum(row) == 0), 'mass_balance_%s' % all_metabolites[i]        
    
        lp.solve()       
        lp.writeLP("flux_mapping.lp")
    
        fluxes_df.loc[all_reactions, PRED_FLUX_L] = \
            map(lambda i: pulp.value(v_pred[i]), all_reactions)
        fluxes_df.loc[measured_reactions, RESID_L] = \
            map(lambda i: pulp.value(v_resid[i]), measured_reactions)
        fluxes_df /= pulp.value(v_pred['Ec_biomass_iJO1366_core_53p95M']), # normalize all fluxes to the biomass flux (i.e. set it to 1)
        
        #%%
        fig, axs = plt.subplots(1, 2, figsize=(14,6))
        fig.subplots_adjust(wspace=0.5)
        axs[0].plot([-50, 50], [-50, 50], 'k', alpha=0.3, linewidth=0.5)
        fluxes_df.plot(kind='scatter', x=MEAS_FLUX_L, y=PRED_FLUX_L, 
                       xerr=MEAS_STDEV_L, ax=axs[0], linewidth=0, s=10,
                       color=(0.7,0.2,0.5))
        axs[0].set_title(condition)
    #    for i in measured_reactions:
    #        xy = fluxes_df.loc[i, [MEAS_FLUX_L, PRED_FLUX_L]]
    #        axs[0].annotate(i, xy, xytext=(10,-5), textcoords='offset points',
    #                        family='sans-serif', fontsize=10, color='darkslategrey')
    
        fluxes_df.loc[~pd.isnull(fluxes_df[RESID_L]), RESID_L].plot(kind='barh', 
                      ax=axs[1], color=(0.7,0.2,0.5))
        axs[1].set_xlabel('residual [mmol/gCDW/h]')
    
        fig.savefig('flux_projection.pdf')
        
        fluxes_df.to_pickle('measured_fluxes.pkl')
        fluxes_df.to_csv('measured_fluxes.csv')
        
        
        #%%
    
    
    #    biomass = model.reactions.get_by_id("Ec_biomass_iJO1366_core_53p95M")
    #    biomass.upper_bound = 0.65
        
        for rid in measured_reactions:
            if fluxes_df["projected fluxes"][rid] < 0 :
                continue
            r = model.reactions.get_by_id(rid)
            r.lower_bound = fluxes_df["projected fluxes"][rid]
            r.upper_bound = fluxes_df["projected fluxes"][rid]   
        
        convert_to_irreversible(model)
        pFBA = optimize_minimal_flux(model, already_irreversible=True)
        reactions = map(str, model.reactions)
        for r in reactions:
            mmol_gCDW_h[condition][r] = pFBA.x_dict[r]

    mmol_gCDW_h = pd.DataFrame.from_dict(mmol_gCDW_h)
    mmol_gCDW_h.to_csv("../data/mmol_gCDW_h.csv")