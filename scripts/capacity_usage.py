# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 15:10:32 2016

@author: dan
"""


from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from load_proteomics import PROTEOMICS
import pandas as pd
import numpy as np
import re

class CAPACITY_USAGE(PROTEOMICS):
    
    def __init__(self):

        self.p = PROTEOMICS()
        self.conditions = self.p.conditions        
        self.cs = self.conditions.index
        self.gr = self.conditions["growth rate [h-1]"]
        self.model = create_cobra_model_from_sbml_file("../data/iJO1366.xml")
        convert_to_irreversible(self.model)
        self.rxns = {r.id:r for r in self.model.reactions}

        mmol_gCDW_h = pd.DataFrame.from_csv("../data/mmol_gCDW_h.csv")
        mmol_gCDW_h = mmol_gCDW_h[self.conditions.index]
        mmol_gCDW_h[mmol_gCDW_h<1e-10] = 0
        self.umol_gCDW_min = mmol_gCDW_h * 1000 / 60 
        self.flux_carrying_rxns = self.flux_carrying_reactions()                
        
        self.enzymes = sorted(map(str, self.model.genes))
        
        self.enzymes_mg_gCDW = self.p.mg_gCDW.copy()
        
        self.enzymes_mg_gCDW.set_index("bnumber", inplace=True)
        self.enzymes_mg_gCDW = self.enzymes_mg_gCDW.loc[self.enzymes, self.cs]        
        self.enzymes_mg_gCDW.replace(np.nan,0, inplace=True)
#        self.enzymes_mg_gCDW.drop("s0001", axis=0, inplace=True)
#        self.enzymes_mg_gCDW = self.enzymes_mg_gCDW.loc[]
        
        
        self.rxns_2_mass = self.reaction_2_mass()
#        


#        self.rxns_per_gene = self.how_many_reactions_per_gene()
        
        self.SA = self.specific_activity()
        self.kmax = self.SA.max(axis=1)
        self.ECU = self.enzyme_capacity_usage()
        self.MC = self.metabolic_capacity()
        self.MCU = self.metabolic_capacity_usage()
                
    def reactions_to_isozymes(self):
        r_to_isozymes = {}
        for rid, r in self.rxns.iteritems():
            isozymes = r.gene_reaction_rule.split("or")
            isozymes = [set(re.findall(r"b[0-9]+", iso)) for iso in isozymes]
            r_to_isozymes[rid] = isozymes
        return r_to_isozymes
        
    def filter_nonexpressed_isozymes(self):
        rxns_to_isozymes = self.reactions_to_isozymes()
        model = create_cobra_model_from_sbml_file("../data/iJO1366.xml")
        reactions = map(str, model.reactions)
        rxns_to_isozymes = {k:v for k, v in rxns_to_isozymes.iteritems()
                            if k in reactions}
        expressed_genes = self.p.copies_fL.copy()
        expressed_genes.set_index("bnumber", inplace=True)

        out = pd.Panel(items=self.cs, 
                       major_axis=self.rxns.keys(), 
                       minor_axis=self.enzymes)
                       
        for c in self.cs:
            # find expressed genes per condition
            # expression threshold is 10 copies per fL
            expression_threshold = 10
            tmp = expressed_genes[c]
            tmp = set(tmp[tmp>expression_threshold].dropna().index) 
            
            for r, iso_list in rxns_to_isozymes.iteritems():

                if r not in self.flux_carrying_rxns[c]:
                    continue
                
                s = set()
                for iso in iso_list:
                    if iso.issubset(tmp):
                        s = s.union(iso)
                out.loc[c, r, list(s)] = 1
#            out = out.to_sparse()
        return out

    def flux_carrying_reactions(self):

        """
            for each condition find the reactions that carry flux
            
            Input: 
                fluxes dataframe
            Output:
                dictionary with conditions as keys and flux carrying reactions
                as values.
        """
        
        out = {}
        for c in self.conditions.index:
            vc = self.umol_gCDW_min[c].copy()
            vc = vc[vc>0]
            out[c] = set(vc.dropna().index) & set(self.rxns.keys())
        return out
    
    def filter_idle_reactions(self, df):
        for c in self.cs:
            for r in df.index:
                if r not in self.flux_carrying_rxns[c]:
                    df.loc[r, c] = np.nan
        return df
                   
    def _enzymes_mg_gCDW(self):
        df = self.p.mg_gCDW.copy()
        df.set_index("bnumber", inplace=True)
        df = df.loc[self.enzymes, self.conditions.index]
        df.replace(np.nan, 0, inplace=True)
        return df
        
        
    def split_non_unique_genes(self, ):
        
        return        
        



                
    
                
    
#        
#            for c in self.cs:
#                if r.id not in self.flux_carrying_rxns[c]:
#                    continue
#                
#            for gene_list in isozymes:
#                for b in gene_list:
#                    if b not in gr_mol.index:
#                        weights.append(subunit_comp[b]*gr_mol[b])
#                    else:
#                        weights.append(np.nan)
#        #    print r.id, weights
#            r_to_weights[r.id] = np.mean(weights)        

    def reaction_2_mass(self):
        out = pd.DataFrame(index=self.rxns.keys(), columns=self.cs)
        out.index.name = "reaction"
        genes = {rid:map(str, r.genes) for rid,r in self.rxns.iteritems()}
        for rid, g_list in genes.iteritems():
            mass = self.enzymes_mg_gCDW.loc[g_list].sum()[self.cs]
            out.loc[rid] = mass
        
        out = self.filter_idle_reactions(out)
        out.replace(0, np.nan, inplace=True)
        return out
        
    def specific_activity(self):
        out = self.umol_gCDW_min.div(self.rxns_2_mass)
        out.replace(0, np.nan, inplace=True)
        return out

    def enzyme_capacity_usage(self):
        return self.SA.div(self.kmax, axis=0)
        
    def metabolic_capacity(self):
        return self.rxns_2_mass.sum()
    
    def metabolic_capacity_usage(self):
        return self.umol_gCDW_min.div(self.kmax, axis=0).sum() / self.MC

if __name__ == '__main__':
    import matplotlib.pyplot as plt    
    cu = CAPACITY_USAGE()
    plt.figure()
    ax = plt.axes()
    plt.scatter(cu.gr, cu.MCU)
    ax.set_xlabel('growth rate [h$^{-1}$]', size=15)
    ax.set_ylabel('capacity usage', size=15)
    [tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]
    ax.set_xlim(0,0.8)
    ax.set_ylim(0,0.8)
    x = cu.filter_nonexpressed_isozymes()
    x.to_pickle("../data/reaction_to_genes_per_condition.pkl")
    y = x.sum()
    y = y[y>1].dropna(how='all')
    mass = cu.enzymes_mg_gCDW.loc[y.index]

    print mass.sum() / cu.enzymes_mg_gCDW.sum()

