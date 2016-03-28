# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 11:57:04 2016

@author: dan
"""

import pandas as pd
import numpy as np

#pd.options.display.precision = 4

class PROTEOMICS(object):
    def __init__(self):

        self.uni_to_b = {row[48:54]:row[0:5].split(';')[0].strip()
                         for row in open("../data/all_ecoli_genes.txt", 'r')}
        self.conditions = pd.DataFrame.from_csv("../data/carbon_sources.csv")
        self.conditions.sort(columns=["growth rate [h-1]"], inplace=True)

        self.copies_cell = pd.read_csv("../data/copies_cell.csv", sep="\t")                             

        self.copies_cell = self._map_upid_2_bnumber(self.copies_cell)
        self.g_mol = self.copies_cell['Molecular weight (Da)'].copy()                
        self.copies_fL = self._copies_cell_2_copies_fL()
        self.mmol_gCDW = self._copies_fL_2_mmol_gCDW()
        self.mg_gCDW = self._mmol_gCDW_2_mg_gCDW()
        self.mass_fraction = self.mass_fraction()
        
    def _map_upid_2_bnumber(self, df):
        
        manual_replacememnts = {
        'D0EX67':'b1107',
        'D4HZR9':'b2755',
        'P00452-2':'b2234',
        'P02919-2':'b0149',
        'Q2A0K9':'b2011',
        'Q5H772':'b1302',
        'Q5H776':'b1298',
        'Q5H777':'b1297',
        'Q6E0U3':'b3183'}
        
        self.uni_to_b.update(manual_replacememnts)
        
        df.replace(to_replace={'UPID':self.uni_to_b}, inplace=True)
        df.drop_duplicates(subset=["UPID"], inplace=True)
        df.rename(columns = {'UPID':'bnumber'}, inplace = True)
        df.sort(columns=["bnumber"], inplace=True) 
        df.index = range(len(df))
        
        return df

    def _copies_cell_2_copies_fL(self):
        df = self.copies_cell.copy()
        volume = self.conditions['single cell volume [fL]']
        copies_cell = df[self.conditions.index]
        copies_fL = copies_cell.div(volume, axis='columns')
        df.loc[df.index, copies_cell.columns] = copies_fL
        
        return df
                
    def _copies_fL_2_mmol_gCDW(self):
        df = self.copies_fL.copy()
        rho = 1100 # average cell density gr/liter
        DW_fraction = 0.3 # fraction of DW of cells
        Avogadro = 6.022 * 1e5 # Avogadro's number
        
        copies_cell = df[self.conditions.index]
        mmol_L = copies_cell / Avogadro
        mmol_gCDW = mmol_L / (rho * DW_fraction)
        df.loc[df.index, copies_cell.columns] = mmol_gCDW
        
        return df
        
    def _mmol_gCDW_2_mg_gCDW(self):

        df = self.mmol_gCDW.copy()
        mmol_gCDW = df[self.conditions.index]        
        mg_gCDW = mmol_gCDW.mul(self.g_mol,axis=0)
        mg_gCDW.replace(np.nan, 0, inplace=True)
        df.loc[df.index, mmol_gCDW.columns] = mg_gCDW
        return df 

    def mass_fraction(self):
        df = self.mg_gCDW.copy()
        mg_gCDW = df[self.conditions.index]
        protein_mass = mg_gCDW.sum()
        mass_fraction = mg_gCDW.div(protein_mass, axis='columns')
        df.loc[df.index, mg_gCDW.columns] = mass_fraction*100
        return df
#    def genes_by_function(self):
#        
#        tree = Tree.FromTMS(open('../data/KO_gene_hierarchy_general.tms', 'r'), 4)
#        f_KEGG = tree.GetNode(name).children
#    
#    
#        reader = csv.reader(open('../data/eco_mapping.csv', 'r'), delimiter='\t')
#        b_to_KEGG = {row[0]:row[2] for row in reader}
#    
#        return {b for b,ko in b_to_KEGG.iteritems() if ko in f_KEGG}

if __name__ == '__main__':
    p = PROTEOMICS()
#    metabolism = p.gene_info[p.gene_info["Annotated functional COG class"] == "METABOLISM"]
#    COG_metabolism = set(metabolism["Annotated functional COG group (description)"])
#    metabolic_mass = p.mg_gCDW.loc[metabolism.index].sum() / p.mg_gCDW.sum() *100

