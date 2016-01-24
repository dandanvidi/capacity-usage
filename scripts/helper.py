import pandas as pd
from trees import Tree
import csv
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
from copy import deepcopy
import numpy as np

protein_info = pd.read_csv('../data/protein_abundance_info.csv', sep='\t')
gc = pd.DataFrame.from_csv("../data/growth_conditions.csv")
gc = gc[gc.reference=='Schmidt et al. 2015']
gr = gc['growth rate [h-1]'][gc.index]
fL_cell = gc['single cell volume [fL]'] /2 # fL (cell volumes are overestimated by a factor of 1.7)
fg_cell_old = pd.read_csv('../data/protein_abundance_[fg_cell].csv')

copies_cell_persist = pd.read_csv('../data/protein_abundance_persistors[copies_cell].csv')

conditions = fL_cell.index & fg_cell_old
new_conditions = copies_cell_persist.columns & fL_cell.index

def map_proteomics(df):
    uni_to_b = {row[48:54]:row[0:5].split(';')[0].strip()
                for row in open("../data/all_ecoli_genes.txt", 'r')}
    
    df.replace(to_replace={'upid':uni_to_b}, inplace=True)
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
    df.replace(to_replace={'upid':manual_replacememnts}, inplace=True)
    df.set_index('upid', inplace=True)                                
    df.index.name = 'bnumber'
    not_identified = ['B8LFD5','D8FH86','D9IX93','E1MTY0','P0CE60','P23477']
    df.drop(not_identified, axis=0, inplace=True)
    df.sort_index(inplace=True)    


def genes_by_function(name):
    
    tree = Tree.FromTMS(open('../data/KO_gene_hierarchy_general.tms', 'r'), 4)
    f_KEGG = tree.GetNode(name).children


    reader = csv.reader(open('../data/eco_mapping.csv', 'r'), delimiter='\t')
    b_to_KEGG = {row[0]:row[2] for row in reader}

    return {b for b,ko in b_to_KEGG.iteritems() if ko in f_KEGG}
    
def convert_copies_fL_to_mmol_gCDW(copies_fL):
 
    rho = 1100 # average cell density gr/liter
    DW_fraction = 0.3 # fraction of DW of cells
    Avogadro = 6.02214129 # Avogadro's number "exponent-less"
    mmol_L = copies_fL / (Avogadro*1e5)
    mmol_gCDW = mmol_L / (rho * DW_fraction)
    return mmol_gCDW
    
def convert_mmol_gCDW_to_mg_gCDW(mmol_gCDW):
        
    protein_info = pd.DataFrame.from_csv('../data/ecoli_genome_info.tsv', sep='\t')
    protein_g_mol = protein_info['molecular_weight[Da]']
    mg_gCDW = mmol_gCDW.mul(protein_g_mol,axis=0)
    mg_gCDW.replace(np.nan, 0, inplace=True)
    return mg_gCDW

def convert_copies_fL_to_mg_gCDW(E):
    tmp = convert_copies_fL_to_mmol_gCDW(E)
    return convert_mmol_gCDW_to_mg_gCDW(tmp)


def get_umol_gCDW_min_from_pFVA(pFVA):

    conds = pFVA.index.levels[0]
    x = pFVA.loc[[(c, 'maximum') for c in conds]]
    x.set_index(conds, inplace=True)
    x = x[x>1e-10]
    return (x * 1000) / 60

def specific_actitivy(V,E,model):

    SA = pd.DataFrame(index=V.index, columns=V.columns)    
    for c in V.columns:
        v = V[c].dropna()
        for rid in v.index:
            if rid == 'flux_counter':
                continue
            r = model.reactions.get_by_id(rid)
            genes = {g.id:g for g in r.genes}
            if 's0001' in genes:
                continue
            weight = []
            for i, (gid, g) in enumerate(genes.iteritems()):
                rxns = {r.id for r in list(g.reactions)} & set(v.index)
                weight.append(E[c][gid] / float(len(rxns)))
            total_weight = sum(weight)
            if total_weight > 0:
                SA[c][rid] = v[rid] / total_weight
    return SA
                                    
def gene_to_flux_carrying_rxns(V,model):
    out = {}
    for c in V.columns:
        out[c] = {}        
        v = V[c].dropna()
        for g in model.genes:
            rxns = {r.id for r in list(g.reactions)} & set(v.index)
            if len(rxns)>0:
                out[c][g.id] = rxns
    return out
    
def get_metabolic_capacity(V,E,model):
    tmp = gene_to_flux_carrying_rxns(V,model)
    capacity = pd.Series({c:E.loc[tmp[c].keys()][c].sum() for c in V.columns})
    return capacity

def get_usage(V,E,model):
    SA = specific_actitivy(V,E,model)
    kmax = SA.max(axis=1)
    return V.div(kmax,axis=0)

def get_efficiency(V,E,model):
    SA = specific_actitivy(V,E,model)
    kmax = SA.max(axis=1)
    return SA.div(kmax,axis=0)
    
def get_capacity_usage(V,E,model):
    capacity = get_metabolic_capacity(V,E,model)
    usage = get_usage(V,E,model)
    return usage.sum() / capacity


def bootstrap_capacity_usage_error(V,E,model,iterations=1000):
    UC = pd.DataFrame(index=range(iterations),columns=V.columns)
    for i in xrange(iterations):
        newE = pd.DataFrame(index=E.index, columns=E.columns)
        for c in V.columns:
            x = E[c]
            x = x[x>0]       
            rand = np.random.choice(x.values, len(x), replace=True)
            newE[c][x.index] = rand
        newE.replace(np.nan, 0, inplace=True)
        UC.loc[i] = get_capacity_usage(V,newE,model)
        print i,
    return UC.std()

def get_foldchange(V,E,gc):
    
    gr = gc['growth rate [h-1]']
    
    combs_all = [(i,j) for (i,j) in combinations(R.gc.index, 2) if gr[j] > gr[i]]    
    delta_mu = pd.Series(data = map(lambda x: np.log2(gr[x[1]]/gr[x[0]]), combs_all),
                         index = combs_all)
    delta_p = pd.DataFrame(index=reactions, columns=combs)
    delta_v = pd.DataFrame(index=reactions, columns=combs)
    for (i, j) in combs:
        delta_p[(i,j)] = np.log2(p[j] / p[i])
        delta_v[(i,j)] = np.log2(v[j] / v[i])
    return delta_p, delta_v, delta_mu

'''            
    x = x.dropna()
    w = w.dropna()
    ix = x.index & w.index
    x = x[ix].values
    w = w[ix].values
    Mw = np.zeros(1000)
    for i in xrange(1000):
        rand = np.random.choice(range(len(x)), len(x), replace=True)
        newx = x[rand]
        neww = w[rand]
        Mw[i] = sum(newx*neww)/sum(neww)
    return np.std(Mw)
'''    
#    print len(fva.keys())    
#    return fva
    
map_proteomics(copies_cell_persist)
map_proteomics(protein_info)
map_proteomics(fg_cell_old)

x = copies_cell_persist[new_conditions]
y = copies_cell_persist['Protein molecular weight']
fg_cell_persist =  x.mul(y,axis=0) / (6.022*1e8) 

fg_cell = fg_cell_old.join(fg_cell_persist, how='outer')
fg_fL = fg_cell.div(fL_cell)

mg_gCDW = fg_fL[gr.index]/(1100/3)*1000 # cell density is 1100 g/L; DW fraction is 1/3    
#mg_gCDW.to_csv('../data/mg_gCDW.csv')
#
out = protein_info.join(mg_gCDW)
#out.to_csv('../data/protein_abundance[mg_gCDW].csv', sep='\t')

#plt.figure()
#ax = plt.axes()
#old = fg_cell_old.index
#new = copies_cell_persist.index
#venn2([old, new], set_labels=('Schmidt et al.', 'Persisters'),set_colors=('#4a6b8a','#801515'),ax=ax)
#plt.tight_layout()
#plt.savefig('../res/comparing coverage.svg')

