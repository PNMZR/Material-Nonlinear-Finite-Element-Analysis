# -*- coding: utf-8 -*-

"""
EX6.6.1.py

Compute element matrix for two-dimensional Laplace equation

Problem description:
    Determine the element matrix for Laplace equation using isoparametric
    four-node quadrilateral element and Gauss-Legendre quadrature for a single
    element.

Variable descriptions:
    k - element matrix
    point2 - integration (or sampling) points
    weight2 - weighting coefficients
    nglx - number of integration points along x-axis
    ngly - number of integration points along y-axis
    xcoord - x coordinate values of nodes
    ycoord - y coordinate values of nodes
    jacob2 - Jacobian matrix
    shape - four-node quadrilateral shape functions
    dhdr - derivatives of shape functions w.r.t. natural coord. r
    dhds - derivatives of shape functions w.r.t. natural coord. s
    dhdx - derivatives of shape functions w.r.t. physical coord. x
    dhdy - derivatives of shape functions w.r.t. physical coord. y
"""

import numpy as np
import feglqd2
import feisoq4
import fejacob2
import federiv2

nnel = 4               # number of nodes per element
ndof = 1               # degrees of freedom per node
edof = nnel * ndof     # degrees of freedom per element

nglx = 2               # use 2x2 integration rule
ngly = 2
xcoord = np.array([-1, 1, 1, -1])              # x coordinate values
ycoord = np.array([-0.75, -0.75, 1.25, 0.25])  # y coordinate values
point2, weight2 = feglqd2.feglqd2(nglx, ngly)  # sampling points & weights

# ---------------------
# numerical integration
# ---------------------
k = np.zeros((edof, edof))            # initialization to zero
for intx in range(nglx):
    x = point2[intx, 0]                # sampling point in x-axis
    wtx = weight2[intx, 0]             # weight in x-axis
    for inty in range(ngly):
        y = point2[inty, 1]            # sampling point in y-axis
        wty = weight2[inty, 1]         # weight in y-axis
        # compute shape functions and derivatives at sampling point
        shape, dhdr ,dhds = feisoq4.feisoq4(x, y)
        # compute Jacobian
        jacob2 = fejacob2.fejacob2(nnel, dhdr, dhds, xcoord, ycoord)
        detjacob = np.linalg.det(jacob2)  # determinant of Jacobian
        invjacob = np.linalg.inv(jacob2)  # inverse of Jacobian matrix
        # derivatives w.r.t.physical coordinate
        dhdx, dhdy = federiv2.federiv2(nnel, dhdr, dhds, invjacob)
        
        # -------------------
        # element matrix loop
        # -------------------
        for i in range(edof):
            for j in range(edof):
                k[i, j] = k[i, j] + (dhdx[i]*dhdx[j] + dhdy[i]*dhdy[j])*wtx*wty*detjacob

print(k)     # print the element matrix