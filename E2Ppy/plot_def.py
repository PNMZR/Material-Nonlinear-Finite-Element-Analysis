# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 16:59:15 2020

@author: lenovo
"""

import numpy as np
import plot_m as plm

def plot_def(fac, u_x, u_y, elemType, dispNodes, dispNodes1, element, node, se):
    """
    Plots the deformed finite element mesh with the support condition
    """
    plm.plot_m(elemType, dispNodes, dispNodes1, element,
               node+fac*np.vstack([u_x, u_y]).T, se)
    #title(' Numerical deformed mesh ')
    # plot(node(dispNodes,1),node(dispNodes,2),'ks')
    # plot(node(dispNodes1,1),node(dispNodes1,2),'ko')