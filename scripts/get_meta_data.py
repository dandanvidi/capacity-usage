# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 16:21:32 2016

@author: dan
"""

from capacity_usage import CAPACITY_USAGE
import pandas as pd
from itertools import product
import numpy as np

gene_info = pd.DataFrame.from_csv('../data/ecoli_genome_info.csv', sep='\t')
copies_fL = pd.DataFrame.from_csv('../data/copies_fL.csv')
flux = pd.DataFrame.from_csv("../data/mmol_gCDW_h.csv")
abundance =  pd.DataFrame.from_csv("../data/g_gCDW.csv")
cu = CAPACITY_USAGE(flux, abundance)


x = pd.DataFrame({r:r.metabolites for r in cu.model.reactions}).T.stack()
r2stoich = x.reset_index()
r2stoich.columns = ['reaction', 'metabolite', 'coefficient']

#%%
l = []
for r,v  in cu.reactions_to_isozymes().iteritems():
    l+=(list(product([r], v, cu.cs)))
r2isozymes = pd.DataFrame(l)
r2isozymes.columns = ['reaction', 'enzyme', 'condition']

#%%
# include abudance of enzyme complexes - take the minimum abudance 
# across all associated genes in the complex
cmplxs = r2isozymes[r2isozymes.enzyme.str.contains(';')].enzyme.drop_duplicates()
for cmplx in cmplxs:
    try: 
        gs = cmplx.split(';')
        copies_fL.loc[cmplx, :] = copies_fL.loc[gs].min(skipna=False)
    except KeyError:
        continue
copies_fL.replace(np.nan, 0, inplace=True)
#%%
#merge reactions and abudances
x = copies_fL.stack()
abundance = x.reset_index()
abundance.columns = ['enzyme', 'condition', 'copies_fL']
df = r2isozymes.merge(abundance, on=['enzyme', 'condition'], how='left')
#
##include also enzyme complexes
#
#for i, j in complexes.iteritems():
#    
#    try:
#        k.loc[k['gene']==j] = copies_fL.loc[gs].min()
#        print gs
#    except KeyError:
#        print i
#%%
#l = [] 
#for r in cu.model.reactions:
#    l+=list(product([r], cu.cs, r.metabolites.keys(), r2isozymes[r.id]))
#df = pd.DataFrame(l, columns=['reaction_ID', 'growth_condition', 'metabolite_ID', 'enzyme'])
#




#df.loc[:, "stoiciometric_coefficient"] = [r.metabolites[m] for r,m in 
#                                          df[['reaction_ID','metabolite_ID']].values]
#df.loc[:, "metabolite_MW[Da]"] = [m.formula_weight for m in df.metabolite_ID]
#
#df.loc[:, "flux[umol_gCDW_min]"] = [cu.umol_gCDW_min.loc[r.id, c] for r,c in 
#                                          df[['reaction_ID','growth_condition']].values]
#
#df.loc[:, "abundance[mg_gCDW]"] = [cu.E.loc[r.id, c] for r,c in 
#                                          df[['reaction_ID','growth_condition']].values]
#                                          
#df.loc[:, "gene_ID"] = [list(r.genes)[0].id if len(r.genes)==1 else None 
#                                for r in df['reaction_ID']]
#
#df.loc[:, "gene_name"] = [gene_info.loc[b, 'gene'] if b in gene_info.index
#                          else None for b in df['gene_ID']]
#
#df.loc[:, "gene_MW[Da]"] = [gene_info.loc[b, 'molecular_weight[Da]'] if b in gene_info.index
#                          else None for b in df['gene_ID']]
#                          
#df.to_csv('../data/df_for_elad.csv', sep='\t')


