# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 13:53:05 2016

@author: dan
"""
import pandas as pd
import re
import numpy as np
from collections import defaultdict
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible

model = create_cobra_model_from_sbml_file("../data/iJO1366.xml")
enzyme_genes = pd.DataFrame.from_csv("../data/model_genes.csv")
gr_mol = enzyme_genes["Molecular weight (Da)"].copy()
convert_to_irreversible(model)

complexes = pd.DataFrame.from_csv("../data/enzyme_complexes.csv")
comp = list(complexes["Gene composition"].values)
comp = [dict(zip(re.findall(r"b[0-9]+", s),re.findall(r"\(([0-9]+)\)", s))) for s in comp]
#
#
all_genes = defaultdict(list)
for s in comp:
    for k,v in s.iteritems():
        all_genes[k].append(float(v))

for b in gr_mol.index:
    if b not in all_genes.keys():
        all_genes[b].append(1.0)
        
subunit_comp = {k:np.mean(v) for k,v in all_genes.iteritems()} 
r_to_weights = {}    
for r in model.reactions:
    isozymes = r.gene_reaction_rule.split("or")
    isozymes = [re.findall(r"b[0-9]+", iso) for iso in isozymes]
    weights = []    
    for genes_list in isozymes:
        for b in genes_list:
            if b in gr_mol.index:
                weights.append(subunit_comp[b]*gr_mol[b])
            else:
                weights.append(np.nan)
#    print r.id, weights
    r_to_weights[r.id] = np.mean(weights)
