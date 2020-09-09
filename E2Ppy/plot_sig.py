# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 09:44:52 2020

@author: lenovo
"""
import matplotlib.pyplot as plt
import plot_stress_field as plsf

def plot_sig(pt1, pt2, pt3, pt4, numx, numy, fac, u_x, u_y, stress, element):
    """
    plots the color coded stress distribution in the finite element
    region along with color bar scale.
    """
  
    fig1, axs1 = plt.subplots(2, 1)
    plsf.plot_stress_field(pt1, pt2, pt3, pt4, numx, numy,
                           fac, u_x, u_y, stress[:,:,0], fig1, axs1[0], element)
    axs1[0].set_title('Stress plot, sigma_xx')
    
    plsf.plot_stress_field(pt1, pt2, pt3, pt4, numx, numy,
                           fac, u_x, u_y, stress[:,:,1], fig1, axs1[1], element)
    axs1[1].set_title('Stress plot, sigma_yy')
    
    fig2, axs2 = plt.subplots(2, 1)
    plsf.plot_stress_field(pt1, pt2, pt3, pt4, numx, numy, fac, 
                    u_x, u_y, stress[:,:,2], fig2, axs2[0], element)
    axs2[0].set_title('Stress plot, sigma_xy')
    
    plsf.plot_stress_field(pt1, pt2, pt3, pt4, numx, numy, fac, 
                    u_x, u_y, stress[:,:,3], fig2, axs2[1], element)
    axs2[1].set_title('Stress plot, sigma_zz')
    

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

