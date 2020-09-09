# -*- coding: utf-8 -*-

import numpy as np

def felp3dt4(x, y, z):
    """
    Purpose:
        Element matrix for three-dimensional Laplace's equation using four-node
        tetrahedral element.

    Synopsis: k = felp3dt4(x, y, z)

    Variable Description:
        k - element matrix (size of 4x4)
        x - x coordinate values of the four nodes
        y - y coordinate values of the four nodes
        z - z coordinate values of the four nodes
    """
    
    k = np.zeros((4, 4))
    xbar = np.array([
        [1, x[0], y[0], z[0]],
        [1, x[1], y[1], z[1]],
        [1, x[2], y[2], z[2]],
        [1, x[3], y[3], z[3]]])
    xinv = np.linalg.inv(xbar)
    vol = (1/6) * np.linalg.det(xbar)

    k[0, 0] = xinv[1, 0]*xinv[1, 0]+xinv[2, 0]*xinv[2, 0]+xinv[3, 0]*xinv[3, 0]
    k[0, 1] = xinv[1, 0]*xinv[1, 1]+xinv[2, 0]*xinv[2, 1]+xinv[3, 0]*xinv[3, 1]
    k[0, 2] = xinv[1, 0]*xinv[1, 2]+xinv[2, 0]*xinv[2, 2]+xinv[3, 0]*xinv[3, 2]
    k[0, 3] = xinv[1, 0]*xinv[1, 3]+xinv[2, 0]*xinv[2, 3]+xinv[3, 0]*xinv[3, 3]
    k[1, 0] = k[0, 1]
    k[1, 1] = xinv[1, 1]*xinv[1, 1]+xinv[2, 1]*xinv[2, 1]+xinv[3, 1]*xinv[3, 1]
    k[1, 2] = xinv[1, 1]*xinv[1, 2]+xinv[2, 1]*xinv[2, 2]+xinv[3, 1]*xinv[3, 2]
    k[1, 3] = xinv[1, 1]*xinv[1, 3]+xinv[2, 1]*xinv[2, 3]+xinv[3, 1]*xinv[3, 3]
    k[2, 0] = k[0, 2]
    k[2, 1] = k[1, 2]
    k[2, 2] = xinv[1, 2]*xinv[1, 2]+xinv[2, 2]*xinv[2, 2]+xinv[3, 2]*xinv[3, 2]
    k[2, 3] = xinv[1, 2]*xinv[1, 3]+xinv[2, 2]*xinv[2, 3]+xinv[3, 2]*xinv[3, 3]
    k[3, 0] = k[0, 3]
    k[3, 1] = k[1, 3] 
    k[3, 2] = k[2, 3] 
    k[3, 3] = xinv[1, 3]*xinv[1, 3]+xinv[2, 3]*xinv[2, 3]+xinv[3, 3]*xinv[3, 3]
    k = vol * k
    return k