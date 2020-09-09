# -*- coding: utf-8 -*-

import numpy as np

def felp2dt3(x, y):
    """
    Purpose:
        Element matrix for two-dimensional Laplace's equation using three-node
        linear triangular element
    
    Synopsis: k = felp2dt3(x, y)

    Variable Description:
        k - element stiffness matrix (size of 3x3)
        x, y - x and y coordinate values of the nodes in the present element
    """
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    y1 = y[0]
    y2 = y[1]
    y3 = y[2]
    k = np.zeros((3, 3))
    # area of the triangle
    A = 0.5 * (x2*y3 + x1*y2 + x3*y1 - x2*y1 - x1*y3 - x3*y2)
    k[0, 0] = ((x3-x2)*(x3-x2) + (y2-y3)*(y2-y3)) / (4*A)
    k[0, 1] = ((x3-x2)*(x1-x3) + (y2-y3)*(y3-y1)) / (4*A)
    k[0, 2] = ((x3-x2)*(x2-x1) + (y2-y3)*(y1-y2)) / (4*A)
    k[1, 0] = k[0, 1]
    k[1, 1] = ((x1-x3)*(x1-x3) + (y3-y1)*(y3-y1)) / (4*A)
    k[1, 2] = ((x1-x3)*(x2-x1) + (y3-y1)*(y1-y2)) / (4*A)
    k[2, 0] = k[0, 2]
    k[2, 1] = k[1, 2]
    k[2, 2] = ((x2-x1)*(x2-x1) + (y1-y2)*(y1-y2)) / (4*A)
    return k