# -*- coding: utf-8 -*-

"""
EX5.10.1.py
To solve the two-dimensional Laplace equation using linear triangular 
elements given as:
    
         u,rr + (u,r)/r + u,zz = 0, 4 < r < 6, 0 < z < 1
         u(4,z) = 100, u,r(6,z) = 20
         u,z(r,0) = O, u,z(r,1) = 0

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
import felpaxt3
import feasmbl1
import feaplyc2

# ---------------------------------
# input data for control parameters
# ---------------------------------
nel = 10               # number of elements
nnel = 3               # number of nodes per element
ndof = 1               # number of dofs per node
nnode = 12             # total number of nodes in system
sdof = nnode * ndof    # total system dofs

# ----------------------------------------------
# input data for nodal coordinate values
# gcoord[i, j] where i-> node No. and j-> r or z
# ----------------------------------------------
gcoord = np.array([
    [4.0, 0.0], [4.0, 1.0], [4.4, 0.0], [4.4, 1.0],
    [4.8, 0.0], [4.8, 1.0], [5.2, 0.0], [5.2, 1.0],
    [5.6, 0.0], [5.6, 1.0], [6.0, 0.0], [6.0, 1.0]])

# ---------------------------------------------------------
# input data for nodal connectivity for each element
# nodes[i, j] where i-> element No. and j-> connected nodes
# ---------------------------------------------------------
nodes = np.array([
    [0,  3, 1], [0,  2,  3],
    [2,  5, 3], [2,  4,  5],
    [4,  7, 5], [4,  6,  7],
    [6,  9, 7], [6,  8,  9],
    [8, 11, 9], [8, 10, 11]])

# ----------------------------------
# input data for boundary conditions
# ----------------------------------
bcdof = np.array([0, 1])        # first and second nodes are constrained
bcval = np.array([100, 100])    # both described value are 100

# --------------------------------------
# initialization of matrices and vectors
# --------------------------------------
ff = np.zeros(sdof)           # initialization of system force vector
kk = np.zeros([sdof, sdof])        # initialization of system matrix
index = np.zeros(nnel*ndof)   # initialization of index vector
ff[10] = 120 * np.pi;              # nodal flux at the outside boundary
ff[11] = 120 * np.pi;              # nodal flux at the outsdie boundary

# --------------------------------------------------------------
# computation of element matrices and vectors and their assembly
# --------------------------------------------------------------
nd = np.zeros(3, dtype='int64')
r = np.zeros(3)
z = np.zeros(3)
for iel in range(nel):
    for i in range(nnel):
        nd[i] = nodes[iel, i]   # connected nodes in (iel)-th element
        r[i] = gcoord[nd[i], 0] # r-coordinates of nodes in (iel)-th element
        z[i] = gcoord[nd[i], 1] # z-coordinates of nodes in (iel)-th element
    # extract system dofs for the element
    index = feeldof.feeldof(nd, nnel, ndof)     
    k = felpaxt3.felpaxt3(r, z)  # compute element matrix
    kk = feasmbl1.feasmbl1(kk, k, index)            # assemble element matrices

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
    r = gcoord[i, 0]
    z = gcoord[i, 1]
    esol[i] = 100-6*20*np.log(4)+6*20*np.log(r)