import sys, os, csv
from trees import Tree
from copy import deepcopy
sys.path.append(os.path.expanduser('~/git/kvivo_max/scripts/'))
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import sys, os
from collections import defaultdict
from catalytic_rates import rates
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from cobra.flux_analysis.single_deletion import single_gene_deletion
from matplotlib_venn import venn3
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation import delete_model_genes
from cobra.core.Gene import parse_gpr, eval_gpr
from collections import Counter
import csv
from collections import defaultdict
sys.path.append(os.path.expanduser('~/git/across-projects'))
from color import ColorMap
from plot_types import cdf
from sklearn import cluster
from uncertainties import unumpy
import textwrap


import pandas as pd
import matplotlib.pyplot as plt

gc = pd.DataFrame.from_csv("../data/growth_conditions.csv")
gc = gc[gc.reference=='Schmidt et al. 2015']
gr = gc['growth rate [h-1]'][gc.index]

data = pd.DataFrame.from_csv('../data/protein_abundance[mg_gCDW].csv', sep='\t')
conditions = data.columns & gr.index

mg_gCDW = data[conditions]

proteins_class = data['Annotated functional COG group (description)']

classes = set(proteins_class.values)

#R = rates()
#
plt.figure(figsize=(8,5))
ax = plt.axes()

b_by_class = {}
tmp = []
for c in classes:
    df = mg_gCDW.loc[proteins_class==c]
    b_by_class[c] = df.index
    if df.sum().max()>20:
        tmp.append(c)

colorlist = ['#8080ff', '#86592d', '#3366ff', '#ffcc00', '#ff33cc', '#00cc66']
colordict = dict(zip(tmp,colorlist))

for c in tmp:
    if type(c)==str:
        label = '\n'.join(textwrap.wrap(c, width=35))
        df = mg_gCDW.loc[proteins_class==c]
        if df.sum().max()>20:
            ax.scatter(gr,df.sum(),s=40,c=colordict[c],label=label,edgecolor='')
            ax.plot(gr, df.sum(), c=colordict[c])
        
ax.set_xlim(0,2)
ax.set_ylim(0)
ax.set_xlabel(r'$\rm{growth\,rate} \left[h^{-1}\right]$', size=fontsize)
ax.set_ylabel(r'$\rm{total\,protein} \left[\frac{\rm{mg}}{\rm{gCDW}}\right]$', size=fontsize)
plt.legend(scatterpoints=1,bbox_to_anchor=(1.65, 1.022))

ax.set_xlabel(r'growth rate $\left[ h^{-1} \right]$', size=fontsize)
ax.set_ylabel(r'total protein $\left[ \frac{mg}{gCDW} \right]$', size=fontsize)

fontsize = 15

[tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]

ax.grid()
#plt.savefig('../res/protein_mass.pdf')
#ufg_cell = unumpy.umatrix(fg_cell, std.values)
#ufg_fl = np.divide(ufg_cell, fL_cell.values)
#umg_gCDW = ufg_fl/(1100/3)*1000
#
#ribo = tree.GetNode('Ribosome')
#ribo = ribo.children
#b = [b for b,ko in b_to_KO.iteritems() if ko in ribo]
#
#ribosome_mass = mg_gCDW.loc[b].dropna(how='all')
#
#plt.figure()
#plt.plot(gr,(ribosome_mass/mg_gCDW.sum()).sum(), 'ro')
#
#
#
#kapp = R.kapp
#kmax = R.kmax['kmax per chain [s^-1]']
#
#efficiency = kapp.div(kmax, axis=0).dropna(how='any')
#efficiency = efficiency[gc.index & efficiency.columns]
#
#
#k=5
#mat = efficiency.as_matrix()
#km = cluster.KMeans(n_clusters=k)
#km.fit(mat)
#labels = km.labels_
#results = pd.DataFrame([efficiency.index,labels]).T
#stack = mat[np.argsort(labels)]
#plt.figure(figsize=(8,10))
#plt.imshow(stack, interpolation='nearest')
#plt.colorbar()
#plt.tight_layout()
#plt.savefig('../res/efficiency.pdf')


#
#fg_cell = pd.read_csv('../data/protein_abundance_[fg_cell].csv', sep='\t')
#fg_cell.replace(to_replace={'upid':uni_to_b}, inplace=True)
#manual_replacememnts = {
#'D0EX67':'b1107',
#'D4HZR9':'b2755',
#'P00452-2':'b2234',
#'P02919-2':'b0149',
#'Q2A0K9':'b2011',
#'Q5H772':'b1302',
#'Q5H776':'b1298',
#'Q5H777':'b1297',
#'Q6E0U3':'b3183'}
#fg_cell.replace(to_replace={'upid':manual_replacememnts}, inplace=True)
#fg_cell.set_index('upid', inplace=True)                                
#fg_cell.index.name = 'bnumber'
#not_identified = ['B8LFD5','D8FH86','D9IX93','E1MTY0','P0CE60','P23477']
#fg_cell.drop(not_identified, axis=0, inplace=True)
#fg_cell.sort_index(inplace=True)    


#SAmax = R.SAmax['max specific activity [umol/mg/min]']
#efficiency = R.SA.div(SAmax,axis='rows')[gc.index]
#
# all expressed enzymes in batch growth in units of mg/gCDW


#expression = pd.DataFrame.from_csv('../data/abundance[copies_fl].csv')
#expression = R._convert_copies_fL_to_mmol_gCDW(expression)
#
#
#def get_functional_group(group_name, extend_fname):
#
#    j = 0    
#    systematic_level = 3
#    genes = []
#    for row in csv.reader(open(extend_fname, 'r'), delimiter='\t'):
#        if len(row) < 3:
#            continue
#        if j != 0:
#            if row[2] != '':
#                break      
#            genes.append(row[systematic_level])        
#  
#        if row[2] == group_name:
#            j = 1
#    return genes
#    
#ribosomal_proteins = get_functional_group('Ribosome', '../data/ecoli_extend.csv')
#
#c = 'GLC_BATCH_mu=0.58_S'
#ribosome_exp = expression.loc[ribosomal_proteins]
#ribosome_exp = ribosome_exp.mean()[gc.index] # mmol/gCDW
#aa_flux = gr * (0.55/110) * 1000 / 3600# mmol/gCDW/s
#f_per_ribo = aa_flux / ribosome_exp # aa/ribosome/s
#
#plt.figure()
#plt.scatter(gr,f_per_ribo)
#plt.xlabel('growth rate [h-1]')
#plt.ylabel('AA per ribosome [s-1]')
#plt.xlim(0)
#plt.ylim(0)
#
#plt.figure()
#plt.plot(gr,ribosome_exp, 'o')
#plt.xlabel('growth rate [h-1]')
#plt.ylabel('ribosomes [mmol/gCDW]')
#plt.xlim(0)
#plt.ylim(0)
#enzymes = R._convert_mmol_gCDW_to_mg_gCDW(R.expression_data[gc.index])
#enzymes.dropna(how='all', inplace=True)
#
#enzymes = enzymes[c].dropna()
#knockouts = set(all_enzymes) - set(enzymes.index)
#knockouts.remove('s0001')
#
##model = create_cobra_model_from_sbml_file('../data/iJO1366.xml')
#model = R.model
#essentiality = single_gene_deletion(model, knockouts)
#essential = set([model.genes.get_by_id(k) for k,v in essentiality[0].iteritems() if v==0])
#
#flux = R.flux_data * 3600
#flux = flux[c]
#tmp = {}
#for g in essential:
#    for r in g.reactions:
#        tmp[g.id] = g.name
#        try:
#            i = flux[r.id] / flux.median()
##            if i > 1e-2:
##                print g
#        except:
#            continue

#
#def support_flux(flux):
#    reactions = [R.rxns[r] for r in flux.index]
#    out = set()
#    for r in reactions:
#        for g in list(r.genes):
#            out.add(g.id)
#    return out
#
#x = enzymes[c].dropna()
#expressed = set(x.index)
#y = flux[c].dropna()
#carry_flux = support_flux(y)
#
#fig = plt.figure()
#ax = plt.axes()
#venn3([all_enzymes, expressed, carry_flux], 
#      ['enzymes\n%s'%len(all_enzymes), 'expressed\n%s'%len(expressed), 'support flux\n%s'%len(carry_flux)], ax=ax)
#
#a = carry_flux - expressed
#out = set()
#i = 0
#for g in a:
#   gene = R.model.genes.get_by_id(g)
#   for r in list(gene.reactions):
#       if 'or' in r.gene_reaction_rule:
#           tree, clean = parse_gpr(r.gene_reaction_rule)
#           i = 1
#           break
#   if i == 1:
#       break
                       
   
 



#util = utilized_enzymes(x,y)
#def utilized_enzymes(expression, flux):
#    out = set()
#    for g in expression.index:
#        e = R.genes[g]
#        reactions = map(lambda x: x.id, e.reactions)
#        if set(reactions) & set(flux.index):
#            out.add(e.id)
#    return out
#    
#for c in gc.index:
#      x = enzymes[c].dropna()
#      expressed = set(x.index)
#      y = flux[c].dropna()
#      carry_flux = support_flux(y)
#      util = utilized_enzymes(x,y)
#      break
#        
#fig = plt.figure()
#ax = plt.axes()
#venn3([all_enzymes, expressed, carry_flux], 
#      ['enzymes\n%s'%len(all_enzymes), 'expressed\n%s'%len(expressed), 'support flux\n%s'%len(util)], ax=ax)
#plt.savefig('../res/enzymes_that_carry_flux.svg')



'''

            

#def perform_pFBA(model, cs, gr, ur):
#
#    rxns = dict([(r.id, r) for r in model.reactions])
#    rxns['EX_glc_e'].lower_bound = 0 # uptake of carbon source reaction is initialized    
#    try:
#        rxns['EX_' + cs + '_e'].lower_bound = -ur # redefine sole carbon source uptake reaction in mmol/gr/h
#    except:
#        print cs, ur
#        rxns['EX_glc_e'].lower_bound = -ur
#    rxns['Ec_biomass_iJO1366_core_53p95M'].upper_bound = gr            
#    print "solving pFBA"
#    solution = optimize_minimal_flux(model, already_irreversible=True)
#    print cs, solution.f
#    flux_dist = pd.DataFrame(model.solution.x_dict.items()).set_index(0)
#    
#    return flux_dist    
#
#fluxes = pd.DataFrame(index=R.rxns.keys(), columns=gc.index)
#for c in gc.iterrows():
#    x = enzymes[c[0]].dropna()
#
#    model = deepcopy(R.model)
#
#    not_expressed = all_enzymes - set(x.index)
#    not_expressed = map(model.genes.get_by_id, not_expressed)
##    delete_model_genes(model, not_expressed)
#
#    cs = c[1]['media_key']
#    gr = c[1]['growth rate [h-1]']
#    ur = c[1]['uptake rate [mmol gCDW-1 h-1]']
#    if np.isnan(ur):
#        ur = 18.5
#    try:        
#        fluxes[c[0]] = perform_pFBA(model, cs, gr, ur)
#    except:
#        print cs
#    break


def bootstrap(x, w):
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

#plt.figure()
#ax = plt.axes()
#for c in gc.index:
#    cdf(enzymes[c], ax=ax)
#    print c, enzymes[c].median()
#ax.set_xscale('log')



unique = set([g for g in R.model.genes if len(g.reactions)==1])
#index = set(proteins.index) & set(map(lambda x:x.id,enzymes))
index = set(proteins.index) & set(map(lambda x:x.id,unique))
expression = proteins.loc[index][gc.index]
mg_gCDW = R._convert_mmol_gCDW_to_mg_gCDW(expression)
mass = pd.DataFrame(index=efficiency.index, columns=gc.index)

for reac in mass.index:
    r = R.rxns[reac]
    genes = map(lambda x: x.id, r.genes)
    try:
        mass.loc[reac] = mg_gCDW.loc[genes].sum()
    except:
        continue
mass.dropna(how='all', inplace=True)

x = gc['growth rate [h-1]'][gc.index]
y = np.zeros(len(gc))
for i,c in enumerate(gc.index):
    a = (mass[c]*efficiency[c]).dropna()
    y[i] = a.sum()/mass.loc[a.index][c].sum()
#
fig = plt.figure(figsize=(8,8))
ax = plt.axes()
(intercept,slope), cov = curve_fit(lambda a,b,x: a*x+b, x, y)
ax.plot(x, slope*x+intercept, 'k:',zorder=0)
for j, i in enumerate(gc.index):
    c = gc.media_key[i]
    mode = gc['growth mode'][i]
    if mode == 'batch':
        ax.scatter(x[j],y[j],c='#ff4d4d',s=80,edgecolor='none',zorder=10)
        ax.errorbar(x[j],y[j],bootstrap(efficiency[i],mass[i]),c='r')
        ax.annotate(c,(x[j],y[j]+0.01),ha='center',va='baseline',size=15)
    else:
        ax.scatter(x[j],y[j], c='y',s=80,edgecolor='none', alpha=0.35)
ax.set_xlabel('growth rate [h$^{-1}$]', size=15)
ax.set_ylabel('effective capacity', size=15)
[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]

plt.grid()
ax.set_xlim(np.floor(10*x.min())/10.,0.72)
ax.set_ylim(np.floor(10*x.min())/10.,0.72)
plt.tight_layout()
plt.savefig('../res/FIG1.png')

fig = plt.figure()
ax = plt.axes()
cm = plt.cm.get_cmap('Blues')
for i,c in enumerate(gc.index):
    a = R.kapp[c].dropna()
#    y = R.SA.loc[a.index][c]
    print a.median()
    cdf(a,color=cm(i/10.),ax=ax,lw=2.5)
kmax = R.kmax['kmax per active site [s-1]'].dropna()
kcat = R.kcat['kcat per active site [s-1]'].dropna()
cdf(kmax,color='k',ax=ax,lw=2.5)
cdf(kcat,color='y',ax=ax,lw=2.5)

ax.set_xscale('log')
ax.set_xlim(1e-2,1e3)

fig = plt.figure()
ax = plt.axes()
ax.plot(x,mass.median(), 'ro')
ax.set_ylim(0,0.2)

#z = {R.rxns[r].id:R.rxns[r].subsystem for r in WCE.index}
#subsystems = defaultdict(list)
#for key, value in sorted(z.iteritems()):
#    subsystems[value].append(key)
#
#colors = ColorMap(subsystems.keys())
#for k,v in subsystems.iteritems():
#    array = matric.loc[v]
#    narray = array.div(array.mean(axis=1), axis=0)
##    narray.dropna(how='any', inplace=True)
#    g = gr[narray.columns]
#    print len(g), len(narray.columns)
#    
#    
##    print k, b
#    ax.plot(g, a*g+b, c=colors[k], marker='o', label=k)
    
#
#fig = plt.figure(figsize=(10,6))
#ax = plt.axes(axisbg='0.95')
#
#plt.scatter(gr, matric.sum(), c='#4DB8FF', edgecolor='none', 
#            s=50)
##            
#labels = [gc['media'][c] for c in gr.index]
#for i, txt in enumerate(labels):
#    ax.annotate(txt, (gr[i],matric.sum()[i]))
##

#plt.grid()

#plt.savefig('../res/growth_rate_and_saturation.pdf')
#x = a[-10:].index
#colors = ColorMap(x)
##plt.figure()
##for r in x:
##    plt.plot(gr, CA.loc[r]/CA.loc[r][0], label=r, c = colors[r], marker='o')
##    plt.legend()
#
#plt.figure()
#for r in x:
#    plt.plot(gr, E.loc[r], label=r, c = colors[r], marker='o')
#    plt.legend()
#    
#plt.figure()
#for r in x:
#    plt.plot(gr, V.loc[r], label=r, c = colors[r], marker='o')
#    plt.legend()
#    plt.figure()
#    plt.hist(a, 40)
#    print a['PSERT'] / gr[c]
#

#plt.tight_layout()
#
##fig = plt.figure(figsize=(6,6))
##ax = plt.axes()
##
##for i, r in enumerate(efficiency.index):
##    x = weighted_concentration.loc[r].astype('float').dropna()
##    y = efficiency.loc[r].astype('float').dropna()
##
##    z = gr[x.index]
##    z.sort
##    x = x[z.index]    
##    y = y[z.index]
##    plt.scatter(x, y, c='b')
##    if i >= 0:
##        break
###ax.set_xscale('log')
##
##
###                
##fig = plt.figure(figsize=(6,6))
##ax = plt.axes()
##conditions = list(efficiency.columns)
###conditions = [conditions[0:12]]# + [conditions[-1]]
##cm = plt.cm.get_cmap('Greens')
###cm = ['r']#, 'b']
##for i, j in enumerate(conditions):
##    x = weighted_concentration[j].astype('float').dropna()
##    y = efficiency[j].astype('float').dropna()
##    
##    index = x.index & y.index
##    x = x.loc[index]
##    y = y.loc[index]
##    
##    plt.scatter(x, y, c=cm(1.0*i/len(conditions)), 
##                edgecolor='none')
##
##[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
##[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]
##ax.set_xlabel('log E [mg/gCDW]', size=15)
##ax.set_ylabel('catalytic efficiency', size=15)
##ax.set_xscale('log')
##ax.set_xlim(1e-4,1e1)
##ax.set_ylim(0,1.1)
##
##
#
##plt.scatter(gr, matric_theoretical, c='#8500AD', edgecolor='none', 
##            s=50, label='relative to $k_{cat}$')
##
##a, p = curve_fit(lambda x, a: a*x, gr, matric_theoretical)
##
###gr = np.append(gr,1/a)
###plt.plot(gr, a* gr)
##
###for i, c in enumerate(R.gc.index):
###    if R.gc['growth mode'][c] == 'batch':
####        plt.scatter(gr[i], matric[c], c='k', edgecolor='none', 
####                    s=50)
###        plt.scatter(gr[i], matric_theoretical[c], c='k', edgecolor='none', 
###                s=50)
##ax.set_xlim(0,1)
##ax.set_ylim(0,0.4)
###plt.legend(scatterpoints=1, loc=3, fontsize=15)#ax.plot([0,0.8],[0,0.8], '#993333')
##plt.grid()
##ax.set_xlabel('growth rate [h$^{-1}$]', size=15)
##ax.set_ylabel('enzyome saturation', size=15)
##
##[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
##[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]
##
##plt.tight_layout()
##
##plt.savefig('../res/mass_efficiency.pdf')
##
##
###fig = plt.figure(figsize=(6,6))
###ax = plt.axes()
###conditions = list(efficiency.columns)
###conditions = [conditions[-10:-1]]# + [conditions[-1]]
###cm = plt.cm.get_cmap('Greens')
###cm = ['b']#, 'b']
###for i, j in enumerate(conditions):
###    x = weighted_concentration[j].astype('float').dropna()
###    y = efficiency[j].astype('float').dropna()
###    
###    index = x.index & y.index
###    x = x.loc[index]
###    y = y.loc[index]
###    
###    plt.scatter(x, y, c=cm[i], 
###                edgecolor='none')
###
###[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
###[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]
###ax.set_xlabel('log E [mg/gCDW]', size=15)
###ax.set_ylabel('catalytic efficiency', size=15)
###ax.set_xscale('log')
###ax.set_xlim(1e-4,1e1)
###ax.set_ylim(0,1.1)
'''