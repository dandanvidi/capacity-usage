# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 21:16:44 2016

@author: dan
"""
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
#from load_proteomics import PROTEOMICS
import pandas as pd
import numpy as np
import re, pickle
import os.path
from collections import Counter
import matplotlib.pyplot as plt

class CAPACITY_USAGE(object):
    
    def __init__(self, flux, abundance):
        
        self.model = create_cobra_model_from_sbml_file("../data/iJO1366.xml")
        convert_to_irreversible(self.model)
        self.rxns = {r.id:r for r in self.model.reactions}
        self._update_subsystems_for_reverse_direction()
        
        g_gCDW = abundance
        mmol_gCDW_h = flux
        mmol_gCDW_h = mmol_gCDW_h[mmol_gCDW_h>1e-10]
        
        self.conditions = pd.DataFrame.from_csv("../data/conditions.csv")        
        self.conditions.sort_values('growth rate Gerosa [h-1]', inplace=True)        
        self.conditions = self.conditions.loc[mmol_gCDW_h.columns &
                                              g_gCDW.columns]        
        self.cs = self.conditions.index
        self.gr = self.conditions["growth rate Gerosa [h-1]"]

        self.cmap = plt.cm.viridis 
        colors = list(self.cmap(np.linspace(0,0.9,6))) + ['#ffd11a']
        self.colors = dict(zip(list(self.cs), colors))         
        
        self.umol_gCDW_min = mmol_gCDW_h[self.cs] * 1e3 / 60
        self.umol_gCDW_min.replace(np.nan, 0, inplace=True)
        self.mg_gCDW = g_gCDW.loc[map(str, self.model.genes), self.cs] * 1e3      
        self.mg_gCDW.replace(np.nan, 0, inplace=True)


        self.flux_carrying_rxns = self.flux_carrying_reactions()                
        self.E = self.map_enzyme_mass_to_reactions()

        self.SA = self.specific_activity()
        self.kmax = self.get_kmax()
        self.vmax = self.get_vmax()
        self.Emin = self.get_Emin()
        self.CU = self.capacity_usage()
        self.MCU = self.metabolic_capacity_usage()
        self.kcat = self.load_kcats()
        
                
    def _update_subsystems_for_reverse_direction(self):
        for rid, r in self.rxns.iteritems():
            if "_reverse" in rid:
                r.subsystem = self.rxns[rid[:-8]].subsystem
                r.name = self.rxns[rid[:-8]].name

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
            vc = self.umol_gCDW_min.loc[:, c]
            out[c] = set(vc[vc>0].index) & set(self.rxns.keys())
        return out

    def reactions_to_isozymes(self):
        r_to_isozymes = {}
        for rid, r in self.rxns.iteritems():
            isozymes = r.gene_reaction_rule.split("or")
            isozymes = [';'.join(re.findall(r"[bs0-9]+", iso)) for iso in isozymes]
            r_to_isozymes[rid] = filter(None, isozymes)
        return r_to_isozymes
        
    def map_reactions_2_enzyme_complexes(self):
        if os.path.isfile("..data/reaction_to_genes_per_condition.pkl"):
            return pickle.load(open("..data/reaction_to_genes_per_condition.pkl","rb"))

        rxns_to_isozymes = self.reactions_to_isozymes()
        out = {}
        for c in self.cs:
            out[c] = pd.DataFrame(index=self.flux_carrying_rxns[c], 
                                  columns=self.mg_gCDW.index)
            # expressed genes in condition "c"
            E_c = self.mg_gCDW.loc[:, c]
            E_c = set(E_c[E_c>0].index)
            
            for r, iso_list in rxns_to_isozymes.iteritems():
                if r not in self.flux_carrying_rxns[c]:
                    continue
                s = set()
                for iso in iso_list:
                    if iso.issubset(E_c):
                        s = s.union(iso)
                    out[c].loc[r, list(s)] = 1
            out[c].dropna(axis=(0,1), how='all', inplace=True)
        pickle.dump(out, open("../data/reaction_to_genes_per_condition.pkl", "wb"))
        return out

    def filter_non_unique_enzymes(self):
        genes = {}
        for c in self.cs:
            genes[c] = []
            for r in self.flux_carrying_rxns[c]:
                g_list = map(str, self.rxns[r].genes)
                genes[c] += g_list
            genes[c] = Counter(genes[c])
            for k, v in genes[c].iteritems():
                if v > 1 :
                    self.mg_gCDW.loc[k, c] = 0 
        return

    def map_enzyme_mass_to_reactions(self):
        '''
            How to distribute enzyme mass between reactions?
            The relation between flux (v) and enzyme amount, or mass (E)
            is not always straightforward as some reactions are catalyzed by
            several enzymes and some enzymes catalyze more than a single reaction.
            We handle this metabolic promescuity in the following manner:
                
            For enzyme i, there exist N reactions it can catalyze and k reactions 
            carry flux (i.e., N>=k and v_j>0 for all j in {1,...,k}
                                       v_l=0 for all l in {k+1,...,N}),
            the enzyme mass associated with reaction j is E/k
            the enzyme mass associated with reaction l is 0
            
            in case k=0 (the enzyme does not carry any reaction, i.e., it is idle):
            the enzyme mass associated with reaction l equals E/N
            
        '''
        if os.path.isfile("../data/reactions_2_mass.csv"):
            return pd.DataFrame.from_csv("../data/reactions_2_mass.csv")
            
        mass = pd.DataFrame(index=self.rxns.keys(), columns=self.cs)
        mass.replace(np.nan, 0, inplace=True)
        for E in self.model.genes:
            reactions = set(map(str, E.reactions))
            for c in self.cs:
                E_mass = self.mg_gCDW.loc[str(E), c]
                inter = reactions.intersection(self.flux_carrying_rxns[c])
                if inter:
                    mass.loc[inter,c] += E_mass / len(inter)
                else:
                    mass.loc[reactions, c] += E_mass / len(reactions)        
        mass.replace(np.nan, 0, inplace=True)
        mass.to_csv("../data/reactions_2_mass.csv")
        return mass
        
    def specific_activity(self):
        out = self.umol_gCDW_min.div(self.E)
        out.replace(np.inf, np.nan, inplace=True)
        return out

    def get_kmax(self):
        return self.SA.max(axis=1)
        
    def load_kcats(self):
        return pd.DataFrame.from_csv("../data/kcat_data.csv")['kcat per AS [s-1]']

    def load_kcats_umolmgmin(self):
        return pd.DataFrame.from_csv("../data/kcat_data.csv")['Vmax [umol min-1 mg-1]']
        
    def kmax_per_sec(self):
        weights = pd.DataFrame.from_csv("../data/kcat_data.csv")['MW per AS [kDa]']
        return self.kmax*weights/60
        
    def get_vmax(self):
        return self.E.mul(self.kmax, axis=0)
        
    def get_Emin(self):
        return self.umol_gCDW_min.div(self.kmax, axis=0)        
        
    def capacity_usage(self):
        vmax = self.vmax.copy()
        vmax.replace(0, 1, inplace=True)
        out = self.umol_gCDW_min.div(vmax, axis=0)
        out.replace(np.inf, np.nan, inplace=True)
        return out 
        
    def metabolic_capacity_usage(self, maximum_rate='kmax'):
        max_rate = self.kmax
        if maximum_rate == 'kcat':
            max_rate = self.kcat
        
        tmp = self.umol_gCDW_min.div(max_rate, axis=0)
        return tmp.sum() / self.mg_gCDW.sum()

    def MCU_stdev(self):
        self.E.mul(self.CU-self.MCU)        
        
        
    def get_master_groups(self):
        master_groups = {'biosynthesis':
                ['Alanine and Aspartate Metabolism',
                 'Arginine and Proline Metabolism',
                 'Cysteine Metabolism',
                 'Glutamate Metabolism',
                 'Glycine and Serine Metabolism',
                 'Histidine Metabolism',
                 'Methionine Metabolism',
                 'Threonine and Lysine Metabolism',
                 'Tyrosine, Tryptophan, and Phenylalanine Metabolism',
                 'Valine, Leucine, and Isoleucine Metabolism',
                 'Murein Biosynthesis',
                 'Lipopolysaccharide Biosynthesis / Recycling',
                 'Cell Envelope Biosynthesis',
                 'Purine and Pyrimidine Biosynthesis',
                 'Cofactor and Prosthetic Group Biosynthesis',
                 ],
                 
                 
                 'central metabolism':
                ['Citric Acid Cycle',
                 'Glycolysis/Gluconeogenesis',
                 'Anaplerotic Reactions',
                 'Pentose Phosphate Pathway',
                 'Oxidative Phosphorylation',
                 'Pyruvate Metabolism'],

                 'rest':
                ['Transport, Inner Membrane',
                 'Transport, Outer Membrane',
                 'Membrane Lipid Metabolism',
                 'Inorganic Ion Transport and Metabolism',
                 'Folate Metabolism',
                 'Unassigned',
                 '',
                 'Glycerophospholipid Metabolism',
                 'Alternate Carbon Metabolism',
                 'Nucleotide Salvage Pathway',
                 'Methylglyoxal Metabolism',
                 'Transport, Outer Membrane Porin',
                 'Murein Recycling',
                 ]}
                 
        return master_groups

def despine(ax, fontsize=15):
    ax.tick_params(right=0, top=0, direction='out', labelsize=fontsize)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel(ax.get_xlabel(), size=15)
    ax.set_ylabel(ax.get_ylabel(), size=15)
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
        
def boot_strap(df, iterations=1000):
    out = np.zeros(len(df.columns))
    for j, c in enumerate(df.columns):
        xc = df[c].dropna().copy()
        median_values = np.zeros(iterations)
        for i in xrange(iterations):
            new_x = np.random.choice(xc, len(xc), replace=True)
            median_values[i] = np.median(new_x)
            out[j] = np.std(median_values)
    return out
    
if __name__ == '__main__':
    flux = pd.DataFrame.from_csv("../data/mmol_gCDW_h.csv")
    abundance =  pd.DataFrame.from_csv("../data/g_gCDW.csv")
    cu = CAPACITY_USAGE(flux, abundance)
    
#%%
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,4))
    cu.MCU.plot(kind='bar', ax=ax1, color='k')
    despine(ax1)
    ax1.set_ylim(0,0.7)
    ax1.set_yticks(np.arange(0,1.1,0.25))
    ax2.scatter(cu.gr, cu.MCU, c='k')
    ax2.set_ylim(0,0.7)
    ax2.set_xlim(0,0.8)
    ax1.set_ylabel('total capacity utilization')
    ax2.set_xlabel('growth rate [$h^{-1}$]')
    ax2.set_xticks(np.arange(0,0.9,0.2))
    ax2.set_yticks(np.arange(0,1.1,0.25))
    despine(ax2)
    plt.tight_layout()
    plt.savefig('../res/MCU_figure.png', dpi=300)
#%%    
    kapp = cu.SA.replace(0, np.nan)
    kapp_m = kapp.median()
    E = cu.E.replace(0, np.nan)
    E_m = E.median()
    
    kapp_ste = boot_strap(kapp)
    E_ste = boot_strap(E)
    
    fig, (ax2, ax1) = plt.subplots(1,2, figsize=(10,4))
    despine(ax1)
    ax1.set_xlim(0,0.8)
    ax1.set_xticks(np.arange(0,0.9,0.2))
    ax1.set_ylim(0,20)
    ax1.set_yticks(np.arange(0,21,4))
    ax1.annotate('A', (0,0.99), xycoords='figure fraction', size=15)
    ax1.set_xlabel(r"growth rate $\left[h^{-1}\right]$")
    ax1.set_ylabel(r"median $k_{app}$ $\left[\frac{\mu mol}{mg \cdot min}\right]$")   
    ax1.errorbar(cu.gr, kapp_m, kapp_ste, fmt='o', color='k')
    

    despine(ax2)
    ax2.set_xlim(0,0.8)
    ax2.set_xticks(np.arange(0,0.9,0.2))
    ax2.set_ylim(0,0.02)
    ax2.set_xlabel(r"growth rate $\left[h^{-1}\right]$")
    ax2.set_ylabel(r"median $E$ $\left[\frac{mg}{gCDW}\right]$") 
    ax2.annotate('B', (0.5,0.99), xycoords='figure fraction', size=15)
    ax2.errorbar(cu.gr, E_m, E_ste, fmt='o', color='k')    
    
    plt.tight_layout()
    
    plt.savefig("kapp and E as a function of growth rate.svg")

