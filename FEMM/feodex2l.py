# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:24:55 2020

@author: lenovo
"""

import numpy as np

def feodex2l(xl, xr, eleng):
    # Purpose:
    # element matrix for (x^2 * u" - 2 * x * u'- 4 * u)
    # using linear element

    # Synopsis:
    # k = feodex2l(xl,xr)

    # Variable Description:
    # k - element matrix (size of 2x2)
    # xl - coordinate value of the left node of the linear element
    # xr - coordinate value of the right node of the linear element

    # element matrix
    
    k = (1/eleng**2)*np.array([[4*xr**2*xl-6*xr*xl**2-xr**3+3*xl**3, 2*xr**2*xl-xr**3-xl**3],
                             [-2*xr*xl**2+xr**3+xl**3, 6*xr**2*xl-4*xr*xl**2-3*xr**3+xl**3]])
    return k