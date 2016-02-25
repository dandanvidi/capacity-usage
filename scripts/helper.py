import cPickle as pickle
import pandas as pd
from trees import Tree
import csv, re
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
from copy import deepcopy
import numpy as np
import seaborn as sb
from collections import defaultdict
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible

ppath = "../../proteomics-collection/"
proteomics = pd.DataFrame.from_csv(ppath+"meta_abundance[copies_fL].csv")
pFBA = pd.DataFrame.from_csv("../data/flux[mmol_gCDW_h].csv")
pFVA = pd.DataFrame.from_csv("../data/flux_variability_[mmol_gCDW_h].csv", header=[0,1]).T
protein_info = pd.read_csv('../data/protein_abundance_info.csv', sep='\t')
gc = pd.DataFrame.from_csv("../data/growth_conditions.csv")
#gc = gc[gc.reference=='Schmidt et al. 2015']
gr = gc['growth rate [h-1]'][gc.index]
fL_cell = gc['single cell volume [fL]'] /2 # fL (cell volumes are overestimated by a factor of 1.7)
fg_cell_old = pd.read_csv('../data/protein_abundance_[fg_cell].csv')

copies_cell_persist = pd.read_csv('../data/protein_abundance_persistors[copies_cell].csv')

model = create_cobra_model_from_sbml_file('../data/iJO1366.xml')
convert_to_irreversible(model)
rxns = {r.id:r for r in model.reactions}

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

def get_complex_molecular_weight(model):
    complexes = pd.DataFrame.from_csv('../data/enzyme_complexes.csv')
    comp = list(complexes['Gene composition'].values)
    comp = [dict(zip(re.findall(r"b[0-9]+", s),re.findall(r"\(([0-9]+)\)", s))) for s in comp]

    protein_info = pd.DataFrame.from_csv('../data/ecoli_genome_info.tsv', sep='\t')
    protein_g_mol = protein_info['molecular_weight[Da]']
    
    all_genes = defaultdict(list)
    for s in comp:
        for k,v in s.iteritems():
            all_genes[k].append(float(v))
    for bnumber in protein_g_mol.index:
        if bnumber not in all_genes.keys():
            all_genes[bnumber].append(1.0)
            
    subunit_comp = {k:np.mean(v) for k,v in all_genes.iteritems()} 
    
    r_to_weights = {}    
    for r in model.reactions:
        isozymes = r.gene_reaction_rule.split('or')
        isozymes = [re.findall(r"b[0-9]+", iso) for iso in isozymes]
        weights = [sum([subunit_comp[b]*protein_g_mol[b] if b in protein_g_mol.index else np.nan 
                    for b in iso]) for iso in isozymes]
        
        r_to_weights[r.id] = np.mean(weights)
    return r_to_weights
    
def convert_copies_fL_to_mg_gCDW(E):
    tmp = convert_copies_fL_to_mmol_gCDW(E)
    return convert_mmol_gCDW_to_mg_gCDW(tmp)


def get_umol_gCDW_min_from_pFVA(pFVA):

    conds = pFVA.index.levels[0]
    x = pFVA.loc[[(c, 'maximum') for c in conds]]
    x.set_index(conds, inplace=True)
    x = x[x>1e-10]
    return (x * 1000) / 60

def gene_to_flux_carrying_rxns(V,model,use_cache=False):
   
    if use_cache:        
        with open('../cache/genes_to_flux_carrying_reactions.p', 'rb') as fp:
            return pickle.load(fp)
    out = {}
    for c in V.columns:
        out[c] = {}        
        vc = V[c]
        vc = vc[vc>0]
        for g in model.genes:
            rxns = {r.id for r in list(g.reactions)} & set(vc.index)
            if len(rxns)>0:
                out[c][g.id] = rxns

    with open('../cache/genes_to_flux_carrying_reactions.p', 'wb') as fp:
        pickle.dump(out, fp)
    return out
    
def convert_SA_to_kcat(SA, MW):
    # MW in units of kDa
    return SA.mul(MW) / 60    
    
    
def flux_carrying_reactions_to_enzymes(V,E,model,use_cache=False):
    if use_cache:       
        with open('../cache/flux_carrying_reactions_to_enzymes.p', 'rb') as fp:
            return pickle.load(fp)
    try:
        V = V.drop('flux_counter')
    except ValueError:
        print "flux couter already removed"
    mapper = {}
    for c in V.columns:
        mapper[c] = {}
        #use only flux carrying reactions in a given condition
        vc = V[c]
        vc = vc[vc>0]
        for rid in vc.index:
            r = model.reactions.get_by_id(rid)
            genes = {g.id:g for g in r.genes}
            # annoing gene in the model - just ignore the reaction it carries
            if 's0001' in genes: continue
            mapper[c][r.id] = {}
            for i, (gid, g) in enumerate(genes.iteritems()):
                rxns = {r.id for r in list(g.reactions)} & set(vc.index)
                mapper[c][rid][gid] = float(len(rxns))
    with open('../cache/flux_carrying_reactions_to_enzymes.p', 'wb') as fp:
        pickle.dump(mapper, fp)
    return mapper

def specific_activity(V,E,model):    
    mapper = flux_carrying_reactions_to_enzymes(V,E,model)
    V = V.to_dict()
    E = E.to_dict()
    SA = {}    
    for c,reactions in V.iteritems():
        SA[c] = {}
        for r,v in reactions.iteritems():
            if r in mapper[c]:
                genes = mapper[c][r]
                abundance = E[c]
                weight = sum([abundance[e] / genes[e] for e in genes])
                if np.isfinite(weight) and weight > 0:
                    SA[c][r] = V[c][r] / weight
                else: 
                    SA[c][r] = np.nan
    SA = pd.DataFrame.from_dict(SA)
    return SA
                                        
def enzyme_capacity_usage(SA):
    kmax = SA.max(axis=1)
    return SA.div(kmax,axis=0)
    
def metabolic_capacity(V,E,model):
    tmp = gene_to_flux_carrying_rxns(V,model)
    capacity = pd.Series({c:E.loc[tmp[c].keys()][c].sum() for c in V.columns})
    return capacity
    
def metabolic_capacity_usage(V,E,model):
    capacity = metabolic_capacity(V,E,model)
    SA = specific_activity(V,E,model)
    ECU = enzyme_capacity_usage(SA)
    E = (V/SA).loc[SA.index]
    return (ECU.mul(E)).sum() / capacity


def bootstrap_capacity_usage_error(V,E,model,iterations=10):
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
    return UC.std()

#def get_foldchange(V,E,gc):
#    
#    gr = gc['growth rate [h-1]']
#    
#    combs_all = [(i,j) for (i,j) in combinations(gc.index, 2) if gr[j] > gr[i]]    
#    delta_mu = pd.Series(data = map(lambda x: np.log2(gr[x[1]]/gr[x[0]]), combs_all),
#                         index = combs_all)
#    delta_p = pd.DataFrame(index=reactions, columns=combs)
#    delta_v = pd.DataFrame(index=reactions, columns=combs)
#    for (i, j) in combs:
#        delta_p[(i,j)] = np.log2(p[j] / p[i])
#        delta_v[(i,j)] = np.log2(v[j] / v[i])
#    return delta_p, delta_v, delta_mu

def get_surface_to_volume_ratio(length,width):
    
    # cylinder + sphere
    volume = np.pi*(length-width)*(width/2)**2 + 4/3*np.pi*(width/2)**3# um^3
    surface = 2*np.pi*(length-width)*(width/2) + 4*np.pi*(width/2)**2# um^2
    
    return surface, volume, surface/volume

def optimize_growth(model, cs):
    
    rxns = {r.id:r for r in model.reactions}
    rxns['EX_glc_e'].lower_bound = 0 # uptake of carbon source reaction is initialized    
    try:
        rxns['EX_' + cs + '_e'].lower_bound = -1000 # redefine sole carbon source uptake reaction in mmol/gr/h
    except KeyError:
        print "%s is not in the model, using glucose instead" %cs
        rxns['EX_glc_e'].lower_bound = -1000
    rxns['Ec_biomass_iJO1366_core_53p95M'].objective_coefficient = 0            
    rxns['Ec_biomass_iJO1366_WT_53p95M'].objective_coefficient = 1            
    model.optimize()
    return   
    
def get_maximal_growth_rate(model, Vmax, condition):
    Vmax = Vmax[condition].copy()
    Vmax = Vmax.dropna()
    Vmax = Vmax * 60 / 1000 # convert to units of mmol/gCDW/h 
    rxns = {r.id:r for r in model.reactions}
    initial_bound = {}
    for r in Vmax.index:
        initial_bound[rxns[r]] = rxns[r].upper_bound
        rxns[r].upper_bound = Vmax[r]
        
    optimize_growth(model, gc['media_key'][condition])
    for r,ub in initial_bound.iteritems():
        r.upper_bound = ub
    return model.solution.f

def get_rand_ECU(ECU,model):
    reactions = [str(r) for r in model.reactions]
    conds = ECU.columns
    rand_ECU = pd.DataFrame(columns=conds, index=reactions)
    for c in conds:
        tmp = ECU[c].dropna()
        rand_ECU[c] = np.random.gamma(tmp.mean(),tmp.std(),len(reactions))
    return rand_ECU

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
    
#map_proteomics(copies_cell_persist)
#map_proteomics(protein_info)
#map_proteomics(fg_cell_old)
#
#x = copies_cell_persist[new_conditions]
#y = copies_cell_persist['Protein molecular weight']
#fg_cell_persist =  x.mul(y,axis=0) / (6.022*1e8) 
#
#fg_cell = fg_cell_old.join(fg_cell_persist, how='outer')
#fg_fL = fg_cell.div(fL_cell)
#
#mg_gCDW = fg_fL[gr.index]/(1100/3)*1000 # cell density is 1100 g/L; DW fraction is 1/3    
##mg_gCDW.to_csv('../data/mg_gCDW.csv')
##
#out = protein_info.join(mg_gCDW)
#out.to_csv('../data/protein_abundance[mg_gCDW].csv', sep='\t')

#plt.figure()
#ax = plt.axes()
#old = fg_cell_old.index
#new = copies_cell_persist.index
#venn2([old, new], set_labels=('Schmidt et al.', 'Persisters'),set_colors=('#4a6b8a','#801515'),ax=ax)
#plt.tight_layout()
#plt.savefig('../res/comparing coverage.svg')

