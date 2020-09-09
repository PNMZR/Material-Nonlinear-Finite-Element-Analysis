# -*- coding: utf-8 -*-

import numpy as np

def federiv2(nnel, dhdr, dhds, invjacob):
    """
    Purpose:
        Determine derivatives of 2-D isoparametric shape functions with respect
        to physical coordinate system.
    
    Synopsis: dhdx, dhdy = federiv2(nnel, dhdr, dhds, invjacob)
    
    Variable Description:
        dhdx - derivative of shape function w.r.t. physical coordinate x
        dhdy - derivative of shape function w.r.t. physical coordinate y
        nnel - number of nodes per element
        dhdr - derivative of shape functions w.r.t. natural coordinate r
        dhds - derivative of shape functions w.r.t. natural coordinate s
        invjacob - inverse of 2-D Jacobian matrix
    """
    dhdx = np.zeros(nnel)
    dhdy = np.zeros(nnel)
    for i in range(nnel):
        dhdx[i] = invjacob[0, 0]*dhdr[i] + invjacob[0, 1]*dhds[i]
        dhdy[i] = invjacob[1, 0]*dhdr[i] + invjacob[1, 1]*dhds[i]
    return dhdx, dhdy