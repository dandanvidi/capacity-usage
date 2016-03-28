# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 12:54:03 2016

@author: dan
"""
import pandas as pd
from capacity_usage import CAPACITY_USAGE
from cobra.io.sbml import create_cobra_model_from_sbml_file

cu = CAPACITY_USAGE()
x = pd.read_pickle("../data/reaction_to_genes_per_condition.pkl")

y = x.sum()
y = y[y>1].dropna(how='all')
mass = cu.enzymes_mg_gCDW.loc[y.index]
model = create_cobra_model_from_sbml_file("../data/iJO1366.xml")
genes = map(cu.model.genes.get_by_id, mass.index)
s = set()
for g in genes:
    s.update(set(map(str,g.reactions)))
print mass.sum() / cu.enzymes_mg_gCDW.sum()