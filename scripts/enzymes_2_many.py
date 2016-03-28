# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 09:49:11 2016

@author: dan
"""

from load_proteomics import PROTEOMICS
from cobra.io.sbml import create_cobra_model_from_sbml_file
from capacity_usage import CAPACITY_USAGE
import pandas as pd

cu = CAPACITY_USAGE()
p = PROTEOMICS()
model = create_cobra_model_from_sbml_file("../data/iJO1366.xml")
enzymes = [b.id for b in model.genes]

flux_carryinf_rxns = cu.flux_carrying_rxns
mg_gCDW = p.mg_gCDW.set_index("bnumber")

enzyme_mass = mg_gCDW.loc[enzymes]

non_unique_mass = pd.DataFrame(index=enzymes, columns=cu.cs)
for c in cu.cs:
    for b in model.genes:
        if set(map(str, b.reactions)).intersection(flux_carryinf_rxns[c]):
            non_unique_mass.loc[b.id, c] = enzyme_mass.loc[b.id, c]

    



#non_unique_mass = mg_gCDW.loc[non_unique]
#
#print non_unique_mass.sum() / enzyme_mass.sum()
