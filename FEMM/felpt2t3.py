# -*- coding: utf-8 -*-

import numpy as np

def felpt2t3(x, y):
    """
    Purpose:
        Element matrix for transient term of two-dimensional Laplace's equation
        using linear triangular element

    Synopsis: m = felpt2t3(x, y)

    Variable Description:
        m - element stiffness matrix (size of 3x3)
        x, y - x and y coordinate values of the nodes in the present element
    """
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    y1 = y[0]
    y2 = y[1]
    y3 = y[2]
    A = 0.5 * (x2*y3 + x1*y2 + x3*y1 -x2*y1 -x1*y3 - x3*y2)
    m = (A/12) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]])
    return m