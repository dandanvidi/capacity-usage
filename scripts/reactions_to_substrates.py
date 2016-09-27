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
flux = pd.DataFrame.from_csv("../data/mmol_gCDW_h.csv")
abundance =  pd.DataFrame.from_csv("../data/g_gCDW.csv")
cu = CAPACITY_USAGE(flux, abundance)
l = [] 
for r in cu.model.reactions:
    l+=list(product([r], cu.cs, r.metabolites.keys()))
df = pd.DataFrame(l, columns=['reaction_ID', 'growth_condition', 'metabolite_ID'])

df.loc[:, "stoiciometric_coefficient"] = [r.metabolites[m] for r,m in 
                                          df[['reaction_ID','metabolite_ID']].values]
df.loc[:, "metabolite_MW[Da]"] = [m.formula_weight for m in df.metabolite_ID]

df.loc[:, "flux[umol_gCDW_min]"] = [cu.umol_gCDW_min.loc[r.id, c] for r,c in 
                                          df[['reaction_ID','growth_condition']].values]

df.loc[:, "abundance[mg_gCDW]"] = [cu.E.loc[r.id, c] for r,c in 
                                          df[['reaction_ID','growth_condition']].values]
                                          
df.loc[:, "gene_ID"] = [list(r.genes)[0].id if len(r.genes)==1 else None 
                                for r in df['reaction_ID']]

df.loc[:, "gene_name"] = [gene_info.loc[b, 'gene'] if b in gene_info.index
                          else None for b in df['gene_ID']]

df.loc[:, "gene_MW[Da]"] = [gene_info.loc[b, 'molecular_weight[Da]'] if b in gene_info.index
                          else None for b in df['gene_ID']]
                          
df.to_csv('../data/df_for_elad.csv', sep='\t')