# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 09:50:08 2020

@author: lenovo
"""
import matplotlib.tri as mtri
import mesh_region as mesh
import matplotlib as mpl
import numpy as np

def plot_field(pt1, pt2, pt3, pt4, numx, numy, fac, u_x, u_y, field, 
               fig, ax, elemType='T3'):
    """
    Forms a color dependent finite element mesh for plotting outputs
    """ 
    node, triangles = mesh.mesh_region(pt1, pt2, pt3, pt4, numx, numy, elemType)
    node = node + fac*np.vstack([u_x, u_y]).T
    
    x = node[:, 0]

    y = node[:, 1]

    z = field

    triang = mtri.Triangulation(x, y, triangles-1)

    # 绘制有限元网格以及等值线图
    p = ax.tricontourf(triang, z, 150, cmap=mpl.cm.jet)
    #ax.triplot(triang, 'ko-')
    fig.colorbar(p, ax=ax)
    ax.axis('equal')
