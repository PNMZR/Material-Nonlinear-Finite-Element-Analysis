# -*- coding: utf-8 -*-

"""
EX5.12.1.py

To solve the three-dimensional Laplace equation for a pyramid shape of domain
using four-node tetrahedral elements. Bottom face has essential boundary
condition and the side faces are insulated.

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
import felp3dt4
import feasmbl1
import feaplyc2

# ---------------------------------
# input data for control parameters
# ---------------------------------
nel = 4                # number of elements
nnel = 4               # number of nodes per element
ndof = 1               # number of dofs per node
nnode = 6              # total number of nodes in system
sdof = nnode * ndof    # total system dofs

# ----------------------------------------------
# input data for nodal coordinate values
# gcoord[i, j] where i-> node No. and j-> x or y
# ----------------------------------------------
gcoord = np.array([
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.5, 0.5],
    [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.5, 0.5, 1.0]])

# ---------------------------------------------------------
# input data for nodal connectivity for each element
# nodes[i, j] where i-> element No. and j-> connected nodes
# ---------------------------------------------------------
nodes = np.array([[3, 0, 2, 5], [0, 1, 2, 5], [1, 4, 2, 5], [4, 3, 2, 5]])

# ----------------------------------
# input data for boundary conditions
# ----------------------------------
bcdof = np.array([0, 1, 3, 4])
bcval = np.array([0, 20, 50, 100])

# --------------------------------------
# initialization of matrices and vectors
# --------------------------------------
ff = np.zeros(sdof)           # initialization of system force vector
kk = np.zeros([sdof, sdof])   # initialization of system matrix
index = np.zeros(nnel*ndof)   # initialization of index vector

# --------------------------------------------------------------
# computation of element matrices and vectors and their assembly
# --------------------------------------------------------------
nd = np.zeros(4, dtype='int64')
x = np.zeros(4)
y = np.zeros(4)
z = np.zeros(4)
for iel in range(nel):
    for i in range(nnel):
        nd[i] = nodes[iel, i]   # connected nodes in (iel)-th element
        x[i] = gcoord[nd[i], 0] # x-coordinates of nodes in (iel)-th element
        y[i] = gcoord[nd[i], 1] # y-coordinates of nodes in (iel)-th element
        z[i] = gcoord[nd[i], 2] # y-coordinates of nodes in (iel)-th element
    # extract system dofs for the element
    index = feeldof.feeldof(nd, nnel, ndof)     
    k = felp3dt4.felp3dt4(x, y, z)          # compute element matrix
    kk = feasmbl1.feasmbl1(kk, k, index)    # assemble element matrices

# -------------------------
# apply boundary conditions
# -------------------------
kk, ff = feaplyc2.feaplyc2(kk, ff, bcdof, bcval)

# -------------------------
# solve the matrix equation
# -------------------------
fsol = np.linalg.solve(kk, ff)