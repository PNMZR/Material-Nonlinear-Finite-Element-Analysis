# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 09:50:08 2020

@author: lenovo
"""
import matplotlib.tri as mtri
import matplotlib as mpl
import numpy as np
import mesh_region as mesh

def plot_stress_field(pt1, pt2, pt3, pt4, numx, numy, fac, u_x, u_y, field, 
               fig, ax, element, elemType='T3'):
    """
    Forms a color dependent finite element mesh for plotting outputs
    """ 
    node, triangles = mesh.mesh_region(pt1, pt2, pt3, pt4, numx, numy, elemType)
    node = node + fac*np.vstack([u_x, u_y]).T
    x = node[:, 0]
    y = node[:, 1]
    triang = mtri.Triangulation(x, y, triangles-1)
    
    z = np.zeros(node.shape[0])
    for n in range(1, node.shape[0]):
        z[n-1] = field[np.where(element == n)].mean()
    # 绘制有限元网格以及等值线图
    p = ax.tricontourf(triang, z, 150, cmap=mpl.cm.jet)
    #ax.triplot(triang, 'ko-')
    fig.colorbar(p, ax=ax)
    ax.axis('equal')
