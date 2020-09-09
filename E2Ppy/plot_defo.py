# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 09:44:52 2020

@author: lenovo
"""
import matplotlib.pyplot as plt
import plot_field as plf

def plot_defo(pt1, pt2, pt3, pt4, numx, numy, fac, 
              u_x, u_y):
    """
    plots the color coded displacement intensity in the finite element 
    region along with color bar scale.
    """
  
    fig, axs = plt.subplots(2, 1)
    field1 = u_x
    plf.plot_field(pt1, pt2, pt3, pt4, numx, numy, fac, 
                   u_x, u_y, field1, fig, axs[0])
    axs[0].set_title('Deformation plot, U_X')
    
    field2 = u_y
    plf.plot_field(pt1, pt2, pt3, pt4, numx, numy, fac, 
                   u_x, u_y, field2, fig, axs[1])
    axs[1].set_title('Deformation plot, U_Y')
    

# function plot_defo(fac,u_x,u_y,elemType)

# % plots the color coded displacement intensity in the finite element
# % region along with color bar scale.

# global node element

# figure
# clf
# subplot(2,1,1);
# plot_field(node+fac*[u_x u_y],element,elemType,u_x);
# colorbar
# title('Deformation plot, U_X')

# subplot(2,1,2);
# plot_field(node+fac*[u_x u_y],element,elemType,u_y);
# colorbar
# title('Deformation plot, U_Y')
# end  % end of function

