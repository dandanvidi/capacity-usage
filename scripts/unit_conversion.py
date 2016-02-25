import numpy as np
import pandas as pd

class Convert(object):

    def __init__(self):
    
        protein_info = pd.DataFrame.from_csv('../data/ecoli_genome_info.tsv', sep='\t')
        self.g_mol = protein_info['molecular_weight[Da]']

    def copies_fL_to_mmol_gCDW(self,copies_fL):
     
        rho = 1100 # average cell density gr/liter
        DW_fraction = 0.3 # fraction of DW of cells
        Avogadro = 6.02214129 # Avogadro's number "exponent-less"
        mmol_L = copies_fL / (Avogadro*1e5)
        mmol_gCDW = mmol_L / (rho * DW_fraction)
        return mmol_gCDW
        
    def mmol_gCDW_to_mg_gCDW(self,mmol_gCDW):
            
        mg_gCDW = mmol_gCDW.mul(self.g_mol,axis=0)
        mg_gCDW.replace(np.nan, 0, inplace=True)
        return mg_gCDW

    def copies_fL_to_mg_gCDW(self,copies_fL):

        tmp = self.copies_fL_to_mmol_gCDW(copies_fL)
        return self.mmol_gCDW_to_mg_gCDW(tmp)

if __name__ == '__main__':
    
    con = Convert()
    ppath = "../../proteomics-collection/"
    copies_fL = pd.DataFrame.from_csv(ppath+"meta_abundance[copies_fL].csv")
    mg_gCDW = con.copies_fL_to_mg_gCDW(copies_fL)
    print mg_gCDW.sum()
