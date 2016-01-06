import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.path.expanduser('~/git/kvivo_max/scripts/'))
from catalytic_rates import rates
from cobra.io.sbml import create_cobra_model_from_sbml_file
import pandas as pd
import numpy as np
from helper import *
from collections import defaultdict

def carry_flux(E, v):
    carry_no_flux = set()
    carry_nonunique_flux = set()
    carry_unique_flux = set()
    E = E[E>0]

    support_flux = set(v.index)
    for i, gid in enumerate(E.index):
        g = genes[gid]
        rxns = {r.id for r in list(g.reactions)}        
        if len(rxns & support_flux) > 1:
            carry_nonunique_flux.add(g.id)
        elif len(rxns & support_flux) == 1:
            carry_unique_flux.add(g.id)
        elif len(rxns & support_flux) == 0:
            carry_no_flux.add(g.id)
            
    return carry_no_flux, carry_unique_flux, carry_nonunique_flux
    
def bnumber_to_functional_group(bnumbers_list=[]):

    reader = csv.reader(open('../data/eco_mapping.csv', 'r'), delimiter='\t')
    b_to_ko = {row[0]:row[2] for row in reader}
    
    ko_gene_hierarchy_fname = '../data/KO_gene_hierarchy_general.tms'    
    systematic_level = 4
    k_to_function = defaultdict(set)    
    for row in csv.reader(open(ko_gene_hierarchy_fname, 'r'), delimiter='\t'):
        if len(row) == systematic_level-1:
            function = row[-1]
        elif len(row) == systematic_level:
            k_to_function[row[-1]].add(function)
    
    out = defaultdict(set)    
    for b in bnumbers_list:
        try:
            function = k_to_function[b_to_ko[b]]
            for f in function:
                out[b].add(f)
        except KeyError:
            continue
            
    return out

R = rates()

model = R.model
genes = {g.id:g for g in model.genes}

ppath = "../../proteomics-collection/"
abundance_copies_fL = pd.DataFrame.from_csv(ppath+"meta_abundance[copies_fL].csv")
abundance_mmol_gCDW = R._convert_copies_fL_to_mmol_gCDW(abundance_copies_fL)

gc = R.gc[R.gc['growth mode']=='batch']
gc = gc[gc.media_key>0]
gc = gc[gc.reference == 'Schmidt et al. 2015']
gc = gc[gc.strain == 'BW25113']

gc.sort(columns=['growth rate [h-1]'], inplace=True)
conditions = gc.index
gr = gc['growth rate [h-1]'].copy()

abundance_mmol_gCDW = abundance_mmol_gCDW[conditions]
protein_info = pd.DataFrame.from_csv('../data/ecoli_genome_info.tsv', sep='\t')
protein_g_mol = protein_info['molecular_weight[Da]']
abundance_mg_gCDW = abundance_mmol_gCDW.mul(protein_g_mol,axis=0)
abundance_mg_gCDW.dropna(how='all', inplace=True)
abundance_mg_gCDW.replace(np.nan, 0, inplace=True)

enzymes_mg_gCDW = abundance_mg_gCDW.loc[genes].dropna(how='all')

enzyme_groups  = pd.DataFrame(index=enzymes_mg_gCDW.index, columns=conditions)
mg_gCDW = pd.DataFrame(index=['no_flux', 'carry_flux', 'unique_flux', 'nonunique_flux'],
                       columns=conditions)

v_mmol_gCDW_s = pd.DataFrame(index=R.rxns.keys(), columns=conditions)
pFVA = pd.DataFrame.from_csv("../data/flux_variability_[mmol_gCDW_s]_relaxation=0.csv", index_col=[0,1])
pVFV = pFVA[pFVA>=0]
pFVA.dropna(how='all', inplace=True)
for i, c in enumerate(conditions):

    vmax = pFVA.loc[c,'maximum']
    v_mmol_gCDW_s[c] = vmax
    may_carry_flux = vmax[vmax>1e-10]

    no_flux, unique_flux, nonunique_flux = carry_flux(enzymes_mg_gCDW[c], may_carry_flux)
    enzyme_groups[c][no_flux] = 0
    enzyme_groups[c][unique_flux] = 1
    enzyme_groups[c][nonunique_flux] = 2

    mg_gCDW[c]['no_flux'] = enzymes_mg_gCDW[c][no_flux].sum()
    mg_gCDW[c]['unique_flux'] = enzymes_mg_gCDW[c][unique_flux].sum()
    mg_gCDW[c]['nonunique_flux'] = enzymes_mg_gCDW[c][nonunique_flux].sum()
    mg_gCDW[c]['carry_flux'] = mg_gCDW[c]['unique_flux'] + mg_gCDW[c]['nonunique_flux']
    
relative_amounts = mg_gCDW.div(enzymes_mg_gCDW.sum())
fig = plt.figure(figsize=(8,5))
fontsize=15
ax = plt.axes()
color = {'no_flux':'lightcoral','carry_flux':'black',
         'unique_flux':'gold','nonunique_flux':'lightskyblue'}
for group in relative_amounts.index:
#    if group in ['no_flux', 'carry_flux']:
    ax.scatter(gr, relative_amounts.loc[group], c=color[group], edgecolor='', label=group)        

ax.set_xlabel(r'growth rate [h$^{-1}$]', size=fontsize)
ax.set_ylabel(r'total mass [%]', size=fontsize)
ax.set_ylim(0,1)
[tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]
plt.legend(scatterpoints=1)
plt.tight_layout()



#for i, c in enumerate(conditions):
#    pFVA = pd.DataFrame.from_csv("../data/flux_variability_%s_[mmol_gCDW_s]_relaxation=0.csv"%c)
##    pFVA[pFVA.maximum<0] == pFVA[pFVA.maximum<0]*(-1) # transport reacions carry negative flux
#    may_carry_flux = pFVA[pFVA.maximum>1e-10]
#    v = v_mmol_gCDW_s[c]
#    v = v[v>0]
#    v.dropna(inplace=True)
#    no_flux, unique_flux, nonunique_flux = carry_flux(enzymes_mg_gCDW[c], v)
#    enzyme_groups[c][no_flux] = 0
#    enzyme_groups[c][unique_flux] = 1
#    enzyme_groups[c][nonunique_flux] = 2
#
#    mg_gCDW[c]['no_flux'] = enzymes_mg_gCDW[c][no_flux].sum()
#    mg_gCDW[c]['unique_flux'] = enzymes_mg_gCDW[c][unique_flux].sum()
#    mg_gCDW[c]['nonunique_flux'] = enzymes_mg_gCDW[c][nonunique_flux].sum()
#    mg_gCDW[c]['carry_flux'] = mg_gCDW[c]['unique_flux'] + mg_gCDW[c]['nonunique_flux']
#    
#relative_amounts = mg_gCDW.div(enzymes_mg_gCDW.sum())
#fig = plt.figure(figsize=(8,5))
#fontsize=15
#ax = plt.axes()
#color = {'no_flux':'lightcoral','carry_flux':'black',
#         'unique_flux':'gold','nonunique_flux':'lightskyblue'}
#for group in relative_amounts.index:
##    if group in ['no_flux', 'carry_flux']:
#    ax.scatter(gr, relative_amounts.loc[group], c=color[group], edgecolor='', label=group)        
#
#ax.set_xlabel(r'growth rate [h$^{-1}$]', size=fontsize)
#ax.set_ylabel(r'total mass [%]', size=fontsize)
#ax.set_ylim(0,1)
#[tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
#[tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]
#plt.legend(scatterpoints=1)
#plt.tight_layout()



no_flux = enzyme_groups[enzyme_groups==0].dropna(how='all').index
b_to_function = bnumber_to_functional_group(no_flux)  
no_flux = pd.DataFrame.from_dict(b_to_function.items()).set_index(0)
#no_flux.to_csv('../res/neveer_carry_flux.tsv', sep='\t')

#    ax.plot(gr, relative_amounts.loc[group], c=color[group])    
#    
#    carry_no_flux_function = protein_info['function'][carry_no_flux]
#    
#    [c] 
#    x = noflux_mg_gCDW / enzymes_mg_gCDW[c].sum()
#    y = uniqueflux_mg_gCDW / enzymes_mg_gCDW[c].sum()
#    z = nonuniqueflux_mg_gCDW / enzymes_mg_gCDW[c].sum()

#
#    ax.plot(gr, y, c='gold')    
#    ax.scatter(gr, y, c='gold', edgecolor='')    
#
#    ax.plot(gr, z, c='lightskyblue')    
#    ax.scatter(gr, z, c='lightskyblue', edgecolor='')    
#    
#

#plt.tight_layout()
#plt.savefig('../res/enzymoe_unique_fraction.pdf')
#
#

#transporters = set()
#tRNA = set()
#for b,v in b_to_function.iteritems():
#    for f in v:
#        if 'transport' in f:
#            transporters.add(b)
#        elif 'tRNA' in f:
#            tRNA.add(b)
            
#new_explained = set([x for x in not_explained 
#                            if b_to_function[x] in 
#                            ['tRNA loading', 
#                            'Insulin signaling pathway', 
#                            'Pertussis', 
#                            'Chaperones and folding catalysts',
#                            'Two-component system',
#                            'Phosphatidylinositol signaling system',
#                            'Sulfur relay system']])    

#SA = pd.DataFrame(index=v_mmol_gCDW_s.index, columns=conditions)
#for i, rid in enumerate(v_mmol_gCDW_s.index):
#    r = R.rxns[rid]
#    genes = map(lambda x: x.id, r.genes)
#    if set(genes).issubset(carry_unique_flux):
#        SA.loc[r.id] = v_mmol_gCDW_s.loc[rid] / enzymes_mg_gCDW.loc[genes].sum()
#
#SA.replace([0, np.inf, -np.inf], np.nan, inplace=True)
#SA.dropna(how='all', inplace=True)
#SA *= (1000 * 60) # convert from mmol/mg/s to umol/mg/min - specific activity
##
##
#mass = np.zeros(len(SA.index))
#g = set()
#for i, rid in enumerate(SA.index):
#    r = R.rxns[rid]
#    genes = map(lambda x: x.id, r.genes)
#    for gid in genes:    
#        g.add(gid)
#    mass[i] = enzymes_mg_gCDW.loc[genes].sum()[c]
#print sum(mass)
#or i, gid in enumerate(carry_unique_flux):
#    g = genes[gid]
    

#for g in carry_unique_flux:
#    genes = map(lambda x: x.id, r.genes)
#    try:
#        SA.loc[r.id] = self.v.loc[r.id] / weighted_mass.loc[genes].sum()
#    except KeyError:
#        continue
#


#return 


#print enzymes_mg_gCDW.loc[carry_unique_flux].sum() / enzymes_mg_gCDW.sum()
#print enzymes_mg_gCDW.loc[carry_nonunique_flux].sum() / enzymes_mg_gCDW.sum()
#def support_flux(flux):
#    reactions = [R.rxns[r] for r in flux.index]
#    out = set()
#    for r in reactions:
#        for g in list(r.genes):
#            out.add(g.id)
#    return out


#reactions = map(lambda x: x.id, R.enzymatic_reactions)
#SA = pd.DataFrame(index=reactions, columns=conditions)
#
#for r in R.enzymatic_reactions.iterkeys():
#    genes = map(lambda g: g.id, r.genes)
#    try:
#        SA.loc[r.id] = v_mmol_gCDW_s.loc[r.id] / enzymes_mg_gCDW.loc[genes].sum()
#    except KeyError:
#        print r
#        continue
#    
#SA.replace([0, np.inf, -np.inf], np.nan, inplace=True)
#SA.dropna(how='all', inplace=True)
#SA *= (1000 * 60)