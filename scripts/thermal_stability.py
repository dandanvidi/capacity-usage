import pandas as pd
from capacity_usage import CAPACITY_USAGE
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import seaborn as sns

#flux = pd.DataFrame.from_csv("../data/mmol_gCDW_h.csv")
abundance =  pd.DataFrame.from_csv("../data/g_gCDW.csv")
#cu = CAPACITY_USAGE(flux, abundance)


uni_to_b = {row[48:54]:row[0:5].split(';')[0].strip()
            for row in open("../data/all_ecoli_genes.txt", 'r')}

id_mapper = pd.DataFrame.from_dict(uni_to_b.items())
id_mapper.columns = ["uniprot", "bnumber"]
TS = pd.read_csv("../data/thermoal_stability_ecoli.csv")
df = TS.merge(id_mapper, on=["uniprot"])
#%%
#cu.mg_gCDW.index.name = "bnumber"
#df.set_index("bnumber")

def despine(ax, fontsize=15):
    ax.tick_params(right=0, top=0, direction='out', labelsize=fontsize)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel(ax.get_xlabel(), size=15)
    ax.set_ylabel(ax.get_ylabel(), size=15)
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    
df = df.merge(abundance, left_on="bnumber", right_index=True)
fig = plt.figure(figsize=(6,6))
ax = plt.axes()

x = df['glucose']
y = df['Protein Abundance']
ax.scatter(x,y, edgecolor='', c='#8080ff')
ax.set_xlim(1e-9, 1e-2)

ax.set_xscale('log')
ax.set_yscale('log')

ax.annotate(r'$r^2=%.2f$'%spearmanr(x,y)[0]**2, (0.1,0.9), 
            xycoords='axes fraction', size=15)

ax.set_xlabel('Schmidt et al. [g/gCDW]', size=15)
ax.set_ylabel('Leuenberger et al. [copies/cell?]', size=15)
ax.tick_params(direction='out', labelsize=15)
plt.tight_layout()