# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 14:47:15 2015

@author: dan
"""

def functional_groups(group='Metabolism'):
    tree = Tree.FromTMS(open('../data/KO_gene_hierarchy_general.tms', 'r'), 4)
    KOs = defaultdict(list)
    for c1 in tree.GetNode(group).children:
        try:
            for c2 in tree.GetNode(c1).children:
        except AttributeError:
            Met_KO[b].append(c2)

            for b in tree.GetNode(c2).children:
                Met_KO[b].append(c2)
                
reader = csv.reader(open('../data/eco_mapping.csv', 'r'), delimiter='\t')
b_to_KO = {row[0]:row[2] for row in reader}
metabolic_genes = {k:Met_KO[v] for k,v in b_to_KO.iteritems() if v in Met_KO}

all_enzymes = {g:metabolic_genes[g] if g in metabolic_genes else '' for g in R.genes}
all_enzymes.pop('s0001')
