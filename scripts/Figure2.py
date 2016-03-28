# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 10:42:32 2016

@author: dan
"""

from sklearn.cluster import KMeans
from capacity_usage import CAPACITY_USAGE
import matplotlib.pyplot as plt
import seaborn as sns
cu = CAPACITY_USAGE()

#%%
data = cu.ECU.dropna(how="any")
subsystems = pd.Series(index=data.index, data=[cu.rxns[r].subsystem for r in data.index])

#%%
AA_biosyn = ['Alanine and Aspartate Metabolism',
             'Arginine and Proline Metabolism',
             'Cysteine Metabolism',
             'Glutamate Metabolism',
             'Glycine and Serine Metabolism',
             'Histidine Metabolism',
             'Methionine Metabolism',
             'Threonine and Lysine Metabolism',
             'Tyrosine, Tryptophan, and Phenylalanine Metabolism',
             'Valine, Leucine, and Isoleucine Metabolism']

central_carbon = ['Citric Acid Cycle',
                  'Glycolysis/Gluconeogenesis',
                  'Anaplerotic Reactions',
                  'Pentose Phosphate Pathway']

central_carbon
#%%
cmap = plt.cm.viridis
plt.figure()
sns.clustermap(data, cmap=cmap, col_cluster=False, metric='chebyshev',
               yticklabels=[])
plt.savefig("../res/cluster_matrix.svg")

#%%
X = cu.ECU.dropna(how='any').values
est = KMeans(n_clusters=11)
est.fit(X)
labels = est.labels_

reacs = list(cu.ECU.dropna(how='any').index)
reacs.sort(key=dict(zip(reacs, labels)).get)

clustered = cu.ECU.loc[reacs].copy()
conds = list(cu.conditions.index)
bg = '0.95'
plt.figure(figsize=(10,10))
ax=plt.axes()
sns.heatmap(clustered, cmap=cmap,ax=ax,xticklabels=conds,yticklabels=False)

r2group = pd.DataFrame(index=reacs, columns=['group','class','value'])
r2group['class'] = labels
r2group['group'] = SUBSYSTEMS.loc[reacs].values
r2group['value'] = np.ones(len(reacs))
groups2class = pd.pivot_table(r2group, index='group', values='value', aggfunc='sum', columns='class').replace(np.nan, 0)

plt.show()
plt.close()