# -*- coding: utf-8 -*-

"""
EX6.6.2.py

To solve the two-dimensional Laplace equation using isoparametric four-node
quadrilateral elements given as:
    
         u,xx + u,yy = 0, 0 < x < 5, 0 < y < 10
         u(x, O) = 0, u(x, 10) = lOOsin(pi*x/10)
         u(O, y) = 0, u,x(5, y) = 0

Variable descriptions:
    k = element matrix
    f = element vector
    kk = system matrix
    ff = system vector
    gcoord = coordinate values of each node
    nodes = nodal connectivity of each element
    index = a vector containing system dofs associated with each element
    bcdof = a vector containing dofs associated with boundary conditions
    bcval = a vector containing boundary condition values associated with
            the dofs in bcdof
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
import feeldof
# import felp2dr4
import feasmbl1
import feaplyc2
import feglqd2

# ---------------------------------
# input data for control parameters
# ---------------------------------
nel = 16               # number of elements
nnel = 4               # number of nodes per element
ndof = 1               # number of dofs per node
nnode = 25             # total number of nodes in system
nglx = 2               # use 2x2 integration rule
ngly = 2               # use 2x2 integration rule
sdof = nnode * ndof    # total system dofs
edof = nnel*ndof       # dofs per element

# ----------------------------------------------
# input data for nodal coordinate values
# gcoord[i, j] where i-> node No. and j-> x or y
# ----------------------------------------------
gcoord = np.array([
    [0.0, 0.0], [1.25, 0.0], [2.5, 0.0], [3.75, 0.0], [5.0, 0.0],
    [0.0, 2.5], [1.25, 2.5], [2.5, 2.5], [3.75, 2.5], [5.0, 2.5],
    [0.0, 5.0], [1.25, 5.0], [2.5, 5.0], [3.75, 5.0], [5.0, 5.0],
    [0.0, 7.5], [1.25, 7.5], [2.5, 7.5], [3.75, 7.5], [5.0, 7.5],
    [0.0, 10.], [1.25, 10.], [2.5, 10.], [3.75, 10.], [5.0, 10.]])

# ---------------------------------------------------------
# input data for nodal connectivity for each element
# nodes[i, j] where i-> element No. and j-> connected nodes
# ---------------------------------------------------------
nodes = np.array([
    [ 0,  1,  6,  5], [ 1,  2,  7,  6], [ 2,  3,  8,  7], [ 3,  4,  9,  8],
    [ 5,  6, 11, 10], [ 6,  7, 12, 11], [ 7,  8, 13, 12], [ 8,  9, 14, 13],
    [10, 11, 16, 15], [11, 12, 17, 16], [12, 13, 18, 17], [13, 14, 19, 18],
    [15, 16, 21, 20], [16, 17, 22, 21], [17, 18, 23, 22], [18, 19, 24, 23]])

# ----------------------------------
# input data for boundary conditions
# ----------------------------------
bcdof = np.array([0, 1, 2, 3, 4, 5, 10, 15, 20, 21, 22, 23, 24])
bcval = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 38.2683, 70.7107, 92.3880, 100])

# --------------------------------------
# initialization of matrices and vectors
# --------------------------------------
ff = np.zeros(sdof)           # initialization of system force vector
kk = np.zeros([sdof, sdof])        # initialization of system matrix
index = np.zeros(nnel*ndof)   # initialization of index vector

# -----------------------------------------------------
# loop for computation and assembly of element matrices
# -----------------------------------------------------
point2, weight2 = feglqd2.feglqd2(nglx, ngly)   # sampling points & weights

nd = np.zeros(4, dtype='int64')
xcoord = np.zeros(4)
ycoord = np.zeros(4)
for iel in range(nel):
    for i in range(nnel):
        nd[i] = nodes[iel, i]   # connected nodes in (iel)-th element
        xcoord[i] = gcoord[nd[i], 0] # x-coordinates of nodes in (iel)-th element
        ycoord[i] = gcoord[nd[i], 1] # y-coordinates of nodes in (iel)-th element
    
    k = np.zeros((edof, edof))  # initialization of element matrix to zero
    
    # ---------------------
    # numerical integration
    # ---------------------
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
            
            # ----------------------
            # compute element matrix
            # ----------------------
            for i in range(edof):
                for j in range(edof):
                    k[i, j] = k[i, j] + (dhdx[i]*dhdx[j] + dhdy[i]*dhdy[j])*wtx*wty*detjacob
    # extract system dofs for the element
    index = feeldof.feeldof(nd, nnel, ndof)   
    # -------------------------
    # assemble element matrices
    # -------------------------
    kk = feasmbl1.feasmbl1(kk, k, index)

# -------------------------
# apply boundary conditions
# -------------------------
kk, ff = feaplyc2.feaplyc2(kk, ff, bcdof, bcval)

# -------------------------
# solve the matrix equation
# -------------------------
fsol = np.linalg.solve(kk, ff)

# -------------------
# analytical solution
# -------------------
esol = np.zeros(nnode)
for i in range(nnode):
    x = gcoord[i, 0]
    y = gcoord[i, 1]
    esol[i] = 100*np.sinh(0.31415927*y)*np.sin(0.31415927*x)/np.sinh(3.1415927)