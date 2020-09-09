# -*- coding: utf-8 -*-

import numpy as np

def felpt2r4(x, y):
    """
    Purpose:
        Element matrix of transient term for two-dimensional Laplace's equation
        using four-node bilinear rectangular element.

    Synopsis: m = felpt2r4(x, y)

    Variable Description:
        m - element stiffness matrix (size of 4x4)
        x, y - x and y coordinate values of the nodes in the present element
    """

    xleng = x[1] - x[0]     # element size in x-axis
    yleng = y[3] - y[0]     # element size in y-axis
    m = (xleng*yleng/36) * np.array([
        [4, 2, 1, 2], [2, 4, 2, 1],
        [1, 2, 4, 2], [2, 1, 2, 4]])
    return m