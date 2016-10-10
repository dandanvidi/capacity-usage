import pandas as pd
from capacity_usage import CAPACITY_USAGE
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import seaborn as sns
from cobra.manipulation.modify import revert_to_reversible
from itertools import product

flux = pd.DataFrame.from_csv("../data/mmol_gCDW_h.csv")
copies_fL = pd.read_csv("../data/abundance[copies_fl].csv")
copies_fL = copies_fL[['bnumber', 'GLC_BATCH_mu=0.58_S']]
abundance = pd.DataFrame.from_csv("../data/g_gCDW.csv")
cu = CAPACITY_USAGE(flux, abundance)


uni_to_b = {row[48:54]:row[0:5].split(';')[0].strip()
            for row in open("../data/all_ecoli_genes.txt", 'r')}

id_mapper = pd.DataFrame.from_dict(uni_to_b.items())
id_mapper.columns = ["uniprot", "bnumber"]
TS = pd.read_csv("../data/thermoal_stability_ecoli.csv")
df = TS.merge(id_mapper, on=["uniprot"])
#%%
model = cu.model.copy()
revert_to_reversible(model)
def one2one_mapping():
    l = []
    for b in cu.model.genes:
        l+=(list(product([b.id], list(b.reactions))))
    df = pd.DataFrame(l)
    df.columns = ['bnumber', 'reaction']
    df.loc[:, '# genes in reaction'] = [len(r.genes) for r in df['reaction']]
    return df

b2r = one2one_mapping()
df = df.merge(b2r, on="bnumber", how='outer')
df = df.merge(copies_fL, on="bnumber", how='outer')
df['reaction'] = df['reaction'].map(str, )
kmax = cu.kmax
kmax.name='kmax'
kcat = cu.load_kcats_umolmgmin()
kcat.name='kcat'
specific_activity = cu.SA.join(kmax)
specific_activity = specific_activity.join(kcat)
specific_activity = specific_activity[['glucose', 'kmax', 'kcat']].dropna(how='all')
specific_activity.columns = ['kapp_glucose', 'kmax', 'kcat']
specific_activity.loc[:, 'kcat/kmax'] = specific_activity['kcat'] / specific_activity['kmax']
df = df.merge(specific_activity, left_on='reaction', right_index=True, how='outer')

df.to_csv("../cache/thermal_stability_metadata.csv")

