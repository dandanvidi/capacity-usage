# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 13:05:46 2016

@author: dan
"""
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
import numpy as np
import seaborn
from capacity_usage import CAPACITY_USAGE


df = pd.DataFrame.from_csv("../data/kcat_data.csv")

MW_AS = df.loc[:, 'MW per AS [kDa]']
kcat  = df.loc[:, 'kcat per AS [s-1]']

fig = plt.figure()
ax = plt.axes()
plt.scatter(kcat, MW_AS**(1/3.))
ax.set_xlim(1e-4, 1e4)
ax.annotate("$R = %.02f$\n$P = %.03f$"%pearsonr(np.log10(kcat), MW_AS**(1/3.)), (1e-3, 6), size=15)
ax.set_xscale('log')
ax.set_xlabel('$k_{cat} [s^{-1}]$', size=20)
ax.set_ylabel('$C \cdot r$', size=20)
ax.set_ylim(1,7)
plt.savefig("../res/peanutbuttermouth.pdf")