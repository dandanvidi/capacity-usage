# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 13:58:55 2016

@author: dan
"""

import os, sys
import matplotlib
from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Wedge, Rectangle

def degree_range(n): 
    start = np.linspace(0,180,n+1, endpoint=True)[0:-1]
    end = np.linspace(0,180,n+1, endpoint=True)[1::]
    mid_points = start + ((end-start)/2.)
    return np.c_[start, end], mid_points
def rot_text(ang): 
    rotation = np.degrees(np.radians(ang) * np.pi / np.pi - np.radians(90))
    return rotation
def gauge(N=5, labels=None, colors='jet_r', 
          cat=1, top_title='', title='', fname='./meter.png'): 
    
    """
    some sanity checks first
    
    """
    
    if not labels:
        labels = ['']*N
    
    if cat > N: 
        raise Exception("\n\nThe category ({}) is greated than the length\nof the labels ({})".format(cat, N))
 
    
    """
    if colors is a string, we assume it's a matplotlib colormap
    and we discretize in N discrete colors 
    """
    
    if isinstance(colors, str):
        cmap = cm.get_cmap(colors)
        cmap = cmap(np.linspace(0,1,N))
        colors = cmap[::-1,:].tolist()
    if isinstance(colors, list): 
        if len(colors) == N:
            colors = colors[::-1]
        else: 
            raise Exception("\n\nnumber of colors {} not equal to number of categories{}\n".format(len(colors), N))

    """
    begins the plotting
    """
    
    fig, ax = plt.subplots()

    ang_range, mid_points = degree_range(N)

    labels = labels[::-1]
    
    """
    plots the sectors and the arcs
    """
    patches = []
    for ang, c in zip(ang_range, colors): 
        # sectors
#        patches.append(Wedge((0.,0.), .4, *ang, facecolor='w', lw=2))
        # arcs
        patches.append(Wedge((0.,0.), .4, *ang, width=0.10, facecolor=c, lw=1, alpha=1))
    
    [ax.add_patch(p) for p in patches]

    """
    set the labels (e.g. 'LOW','MEDIUM',...)
    """

    for mid, lab in zip(mid_points, labels): 

        ax.text(0.35 * np.cos(np.radians(mid)), 0.35 * np.sin(np.radians(mid)), lab, \
            horizontalalignment='center', verticalalignment='center', fontsize=25, \
            fontweight='bold', rotation = rot_text(mid))

    """
    set the bottom banner and the title
    """
    r = Rectangle((-0.45,-0.),0.9,0.001, facecolor='w', lw=2)
    ax.add_patch(r)
#    ax.line()
    ax.text(0, -0.08, title, horizontalalignment='center', \
         verticalalignment='center', fontsize=22, fontweight='bold')

    """
    plots the arrow now
    """
    
    pos = mid_points[np.abs(cat*N - N)]
    
    ax.arrow(0, 0, 0.225 * np.cos(np.radians(pos)), 0.225 * np.sin(np.radians(pos)), \
                 width=0.02, head_width=0.05, head_length=0.1, fc='k', ec='k')
    
#    ax.plot([0, 0.015], [0.225* np.cos(np.radians(pos)),0.225 * np.sin(np.radians(pos))], c='g', lw=2)
#    ax.plot([0, -0.015], [0.225,0], c='k', lw=2)
    ax.add_patch(Circle((0, 0), radius=0.02, facecolor='k'))
    ax.add_patch(Circle((0, 0), radius=0.01, facecolor='w', zorder=11))

    """
    removes frame and ticks, and makes axis equal and tight
    """
    
    ax.set_frame_on(False)
    ax.axes.set_xticks([])
    ax.axes.set_yticks([])
    ax.axis('equal')
    ax.set_title(top_title, fontsize=25)
    plt.tight_layout()

    fig.savefig(fname)

if __name__ == '__main__':
    N = 30
    gauge(N=N, colors='viridis', cat=0.68, 
          title=r'', fname='gauge.svg') 
          
    