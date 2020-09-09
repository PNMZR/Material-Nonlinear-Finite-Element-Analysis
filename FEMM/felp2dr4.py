# -*- coding: utf-8 -*-

import numpy as np

def felp2dr4(x, y):
    """
    Purpose:
        Element matrix for two-dimensional Laplace's equation using 
        four-node bilinear rectangular element

    Synopsis: k = felp2dr4(x, y)

    Variable Description:
        k - element stiffness matrix (size of 4x4)
        x, y - x and y coordinate values of the nodes in the present element
    """
    
    xleng = x[1] - x[0]     # element size in x-axis
    yleng = y[3] - y[0]     # element size in y-axis
    k = np.zeros((4, 4))
    k[0, 0] = (xleng*xleng + yleng*yleng) / (3*xleng*yleng)
    k[0, 1] = (xleng*xleng - 2*yleng*yleng) / (6*xleng*yleng)
    k[0, 2] = -0.5 * k[0, 0]
    k[0, 3] = (yleng*yleng - 2*xleng*xleng) / (6*xleng*yleng)
    k[1, 0] = k[0, 1]
    k[1, 1] = k[0, 0]
    k[1, 2] = k[0, 3]
    k[1, 3] = k[0, 2]
    k[2, 0] = k[0, 2]
    k[2, 1] = k[1, 2]
    k[2, 2] = k[0, 0]
    k[2, 3] = k[0, 1]
    k[3, 0] = k[0, 3]
    k[3, 1] = k[1, 3]
    k[3, 2] = k[2, 3]
    k[3, 3] = k[0, 0]
    return k