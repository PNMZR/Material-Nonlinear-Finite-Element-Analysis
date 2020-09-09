# -*- coding: utf-8 -*-

"""
EX5.9.1.py

To solve the two-dimensional Laplace equation using linear triangular 
elements given as:
    
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
"""

import numpy as np
import feeldof
import felp2dt3
import feasmbl1
import feaplyc2

# ---------------------------------
# input data for control parameters
# ---------------------------------
nel = 32               # number of elements
nnel = 3               # number of nodes per element
ndof = 1               # number of dofs per node
nnode = 25             # total number of nodes in system
sdof = nnode * ndof    # total system dofs

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
    [ 0,  1,  6], [ 1,  2,  7], [ 2,  3,  8], [ 3,  4,  9],
    [ 0,  6,  5], [ 1,  7,  6], [ 2,  8,  7], [ 3,  9,  8],
    [ 5,  6, 11], [ 6,  7, 12], [ 7,  8, 13], [ 8,  9, 14],
    [ 5, 11, 10], [ 6, 12, 11], [ 7, 13, 12], [ 8, 14, 13],
    [10, 11, 16], [11, 12, 17], [12, 13, 18], [13, 14, 19],
    [10, 16, 15], [11, 17, 16], [12, 18, 17], [13, 19, 18],
    [15, 16, 21], [16, 17, 22], [17, 18, 23], [18, 19, 24],
    [15, 21, 20], [16, 22, 21], [17, 23, 22], [18, 24, 23]])

# ----------------------------------
# input data for boundary conditions
# ----------------------------------
bcdof = np.array([0, 1, 2, 3, 4, 5, 10, 15, 20, 21, 22, 23, 24])
bcval = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 38.2683, 70.7107, 92.3880, 100])

# --------------------------------------
# initialization of matrices and vectors
# --------------------------------------
ff = np.zeros([sdof, 1])           # initialization of system force vector
kk = np.zeros([sdof, sdof])        # initialization of system matrix
index = np.zeros([nnel*ndof, 1])   # initialization of index vector

# --------------------------------------------------------------
# computation of element matrices and vectors and their assembly
# --------------------------------------------------------------
nd = np.zeros(3, dtype='int64')
for iel in range(nel):
    nd[0] = nodes[iel, 0]  # 1st connected node for (iel)-th element
    nd[1] = nodes[iel, 1]  # 2nd connected node for (iel)-th element
    nd[2] = nodes[iel, 2]  # 3rd connected node for (iel)-th element
    # coord values of 1st node
    x1 = gcoord[nd[0], 0]
    y1 = gcoord[nd[0], 1]
    # coord values of 2nd node
    x2 = gcoord[nd[1], 0]
    y2 = gcoord[nd[1], 1]
    # coord values of 3rd node
    x3 = gcoord[nd[2], 0]
    y3 = gcoord[nd[2], 1]
    # extract system dofs for the element
    index = feeldof.feeldof(nd, nnel, ndof)     
    k = felp2dt3.felp2dt3(x1, y1, x2, y2, x3, y3) # compute element matrix
    kk = feasmbl1.feasmbl1(kk, k, index)          # assemble element matrices

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