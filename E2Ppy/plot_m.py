# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 23:14:04 2020

@author: lenovo
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_m(elemType, dispNodes, dispNodes1, element, node, se, deformed_mesh=False):
    """
    Plots the finite element mesh with the support condition
    global node element
    """
    fig, ax = plt.subplots()
    if deformed_mesh == False:
        ax.set_title('Undeformed FE mesh')
    else:
        ax.set_title('Deformed FE mesh')
    for e in range(1, element.shape[0]+1):
        if elemType == 'Q9':                   # 9-node quad element
            ord = [1, 5, 2, 6, 3, 7, 4, 8, 1]
        elif elemType == 'Q8':                 # 8-node quad element
            ord = [1, 5, 2, 6, 3, 7, 4, 8, 1]
        elif elemType == 'T3':                 # 3-node triangle element
            ord = [1, 2, 3, 1]
        elif elemType == 'T6':                 # 6-node triangle element
            ord = [1, 4, 2, 5, 3, 6, 1]
        elif elemType == 'Q4':                  # 4-node quadrilateral element
            ord = [1, 2, 3, 4, 1]
        elif elemType == 'L2':                  # 2-node line element
            ord = [1, 2]   
  
        xpt = np.zeros(len(ord))
        ypt = np.zeros(len(ord))
        for n in range(1, len(ord)+1):
            xpt[n-1] = node[element[e-1, ord[n-1]-1]-1, 0]
            ypt[n-1] = node[element[e-1, ord[n-1]-1]-1, 1]     
        ax.plot(xpt, ypt, se)
    ax.plot(node[dispNodes-1, 0], node[dispNodes-1, 1], 'ks')
    ax.plot(node[dispNodes1-1, 0], node[dispNodes1-1, 1], 'ko')
    ax.axis('off')
