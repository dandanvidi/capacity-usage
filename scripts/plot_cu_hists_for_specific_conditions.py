# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 09:03:53 2016

@author: dan
"""

import pandas as pd
from capacity_usage import CAPACITY_USAGE
import matplotlib.pyplot as plt
#import seaborn as sns 

import sys, os
import numpy as np
from scipy.stats import ranksums, wilcoxon
import re

cmap = plt.cm.viridis

flux = pd.DataFrame.from_csv("../data/mmol_gCDW_h.csv")
abundance =  pd.DataFrame.from_csv("../data/g_gCDW.csv")
cu = CAPACITY_USAGE(flux, abundance)


#%%

conds = cu.conditions.index
color = ['red', '0.7']
colors = dict(zip(conds,color))
#hpad=0.5
fs=15
figsize=(8,6)
b = 10
for c in conds:
    capacity_usage = cu.CU.loc[:, c].copy()
    capacity_usage.replace(lambda x: x<0, np.nan, inplace=True)
    capacity_usage.dropna(inplace=True)    
    out, bins = pd.cut(capacity_usage,b, retbins=True)
    out.name = "bins"
    df = pd.DataFrame(out)
    s = df.groupby("bins")
#    E = cu.E.loc[capacity_usage.index, c] / cu.E.loc[capacity_usage.index, c].sum() 
    cap = {k:cu.E.loc[v, c].sum() for k,v in s.groups.iteritems()}
    cap = {k:v/sum(cap.values()) for k,v in cap.iteritems()}
    cap  = pd.Series(cap)

#    df = capacity_usage/2
#    print len(df)
#    hist, bin_edges = np.histogram(df, bins=50, normed=True)
#    bin_edges *= 2
    C = [cmap(1-i) for i in bins]
#    ymax = hist.max() 
#    print hist[-1]
#    ylim=15    
#    if ymax < ylim:
#
    plt.figure(figsize=figsize)
    ax = plt.axes()
#        ax = plt.subplot2grid((4, 2), (0, 0), colspan=2, rowspan=4)
#    ax.bar(bin_edges[:-1],hist,color=C,width=0.02)
    cap.plot(kind='bar', width=0.75, ax=ax, color=C)          
    ax.set_ylim(0,.40)
    ax.tick_params(right=0, top=0, bottom=0, direction='out', labelsize=fs)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_yticks(np.arange(0,.41,.10))

    ylabels = []
    for i in cap.index:
        tmp = re.findall("[-+]?\d*\.\d+|\d+", i)
        ylab = "-".join(tmp)
        if "-0.001" in ylab:
            ylab = '0-0.1'
        ylabels.append(ylab)
    ax.set_xticklabels(ylabels, rotation=45)    
    ax.set_xlabel('capacity utilization', size=fs)
    ax.set_ylabel('fraction of total enzyme mass', size=fs)
    ax.annotate(c, (4, 0.35), size=fs, color='0.5')
    plt.tight_layout()
    plt.savefig('../res/hist_%s.svg'%c)    
#    ax.spines['bottom'].set_visible(False)

'''
        hpad=None
    else:

        plt.figure(figsize=figsize)
        ax_up = plt.subplot2grid((4, 2), (0, 0), colspan=2, rowspan=1)
        ax_low = plt.subplot2grid((4, 2), (1, 0), colspan=2, rowspan=3)
        
        # set upper axis properties
        ax = ax_up

        ax.set_xticklabels([])
        ax.set_ylim(45,60)
        ax.set_xlim(0,1)
        ax.set_yticks([40, 60])
        ax.set_yticklabels(['', 60])
        d = .015
        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        ax.plot((-d,+d),(-d,+d), **kwargs)
        
        ax.bar(bin_edges[:-1],hist,color=C,width=0.02)
        ax.set_ylabel('')
        
        ax = ax_low
        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        ax.plot((-d,+d),(1-d,1+d), **kwargs)
        ax.bar(bin_edges[:-1],hist,color=C,width=0.02)
        ax.set_ylim(0,ylim)
        ax.set_yticks(np.arange(0, ylim-2, 5))
        ax.set_yticklabels([0,5,10,''])
    # set lower axis properties
        
    ax.set_ylabel('')
    ax.tick_params(right=0, top=0, direction='out', labelsize=fs)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xlim(0,1)
    plt.tight_layout(h_pad=hpad)
    
    plt.savefig('../res/hist_%s.svg'%c)
#%%
'''