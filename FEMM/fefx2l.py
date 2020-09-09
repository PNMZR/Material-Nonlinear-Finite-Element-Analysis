# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:24:55 2020

@author: lenovo
"""

import numpy as np

def fefx2l(xl, xr, eleng):
    # Purpose:
    # element vector for f(x)=x^2
    # using linear element

    # Synopsis:
    # f = fefx2l(xl, xr, eleng)

    # Variable Description:
    # f - element vector (size of 2xl)
    # xi - coordinate value of the left node
    # xr - coordinate value of the right node

    # element vector
    k = (1/12*eleng)*np.array([-4*xr*xl**3+xr**4+3*xl**4, -4*xr**3*xl+3*xr**4+xl**4])
    return k