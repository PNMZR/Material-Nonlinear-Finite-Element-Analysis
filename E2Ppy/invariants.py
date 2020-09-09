# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 11:44:35 2020

@author: lenovo
"""
import numpy as np

def invariants(stress):
    """
    calculates the stress invariants for a single gauss point
    """
    s_xx = stress[0]
    s_yy = stress[1]
    t_xy = stress[2]
    s_zz = stress[3]
    p = (s_xx + s_yy + s_zz)/3
    t = np.sqrt((s_xx-s_yy)**2+(s_yy-s_zz)**2+(s_zz-s_xx)**2+6*t_xy**2)/np.sqrt(3)
    q = np.sqrt(1.5)*t
    sx = s_xx - p
    sy = s_yy - p
    sz = s_zz - p

    if q < 1e-6:
       theta = 0
    else:   
       j3 = sx*sy*sz-(sx*t_xy**2)
       sine = -3*j3*np.sqrt(6)/t**3
       if sine > 1:
           sine = 1
       if sine < -1:
           sine = -1
       theta = 1/3*(np.arcsin(sine)/np.pi*180)
    return p, q, theta