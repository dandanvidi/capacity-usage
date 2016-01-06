import pandas as pd
from trees import Tree
import csv
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
from copy import deepcopy


protein_info = pd.read_csv('../data/protein_abundance_info.csv', sep='\t')
gc = pd.DataFrame.from_csv("../data/growth_conditions.csv")
gc = gc[gc.reference=='Schmidt et al. 2015']
gr = gc['growth rate [h-1]'][gc.index]
fL_cell = gc['single cell volume [fL]']/2 # fL (cell volumes are overestimated by a factor of 1.7)
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
    
from cobra.core import Metabolite, Reaction
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from cobra.flux_analysis.variability import flux_variability_analysis

def perform_pFVA(model, cs, gr, ur, reactions, fraction_of_optimum=1.0):

    model = deepcopy(model)
    rxns = dict([(r.id, r) for r in model.reactions])

    rxns['EX_glc_e'].lower_bound = 0 # uptake of carbon source reaction is initialized    
    rxns['EX_' + cs + '_e'].lower_bound = -ur # redefine sole carbon source uptake reaction in mmol/gr/h

    rxns['Ec_biomass_iJO1366_core_53p95M'].upper_bound = gr        
    rxns['Ec_biomass_iJO1366_core_53p95M'].lower_bound = gr     
#    rxns['PDX5PO2'].upper_bound = 0
    fake = Metabolite(id='fake')
    model.add_metabolites(fake)        
            
    for r in model.reactions:
        r.add_metabolites({fake:1})
        
    flux_counter = Reaction(name='flux_counter')
    flux_counter.add_metabolites(metabolites={fake:-1})                

    model.add_reaction(flux_counter) 
    model.change_objective(flux_counter)
    
    fva = flux_variability_analysis(model, 
                                    reaction_list=reactions, 
                                    objective_sense='minimize',
                                    fraction_of_optimum=fraction_of_optimum)
                                    
    return fva
    
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

