# -*- coding: utf-8 -*-

import numpy as np

def fejacob2(nnel, dhdr, dhds, xcoord, ycoord):
    """
    Purpose:
        Determine the Jacobian for two-dimensional mapping.
        
    Synopsis: jacob2 = fejacob2(nnel, dhdr, dhds, xcoord, ycoord)
    
    Variable description:
        jacob2 - Jacobian for one-dimension
        nnel - number of nodes per element
        dhdr - derivative of shape functions w.r.t. natural coordinate r
        dhds - derivative of shape functions w.r.t. natural coordinate s
        xcoord - x axis coordinate values of nodes
        ycoord - y axis coordinate values of nodes
    """

    jacob2 = np.zeros((2, 2))
    for i in range(nnel):
        jacob2[0, 0] = jacob2[0, 0] + dhdr[i]*xcoord[i]
        jacob2[0, 1] = jacob2[0, 1] + dhdr[i]*ycoord[i]
        jacob2[1, 0] = jacob2[1, 0] + dhds[i]*xcoord[i]
        jacob2[1, 1] = jacob2[1, 1] + dhds[i]*ycoord[i]
    return jacob2