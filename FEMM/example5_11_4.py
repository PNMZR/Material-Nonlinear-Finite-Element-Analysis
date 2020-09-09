# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import feeldof
import felp2dr4
import felpt2r4
import feasmbl1
import feaplyc2
import feflxl2
import fef1l
import feasmbl2

"""
EX5.11.4.py
To solve the transient two-dimensional Laplace's equation using bilinear 
rectangular elements and backward difference method, given as:
    
           a*u,t = u,xx + u,yy  0 < x < 0.02, 0 < y < 0. 01
           
Boundary conditions:
    
           u(0,y,t) = 300, u(0.02,y,t) = 300
           u,y(x,0,t) = 0, u,y(x,0.01,t) = 20(u-50)
           
Initial condition:
    
           u(x,y,0) = 300 over the domain

Variable descriptions
    k = element matrix for time-independent term (u,xx + u,yy)
    m = element matnx for time-dependent term (u,t)
    f = element vector
    kk = system matrix of k
    mm = system matrix of m
    ff = system vector
    fn = effective system vector
    fsol = solution vector
    sol = time history solution of selected nodes
    gcoord = coordinate values of each node
    nodes = nodal connectivity of each element
    index = a vector containing system dofs associated with each element
    bcdof = a vector containing dofs associated with boundary conditions
    bcval = a vector containing boundary condition values associated with
            the dofs in bcdof
    k1 = element matrix due to Cauchy-type flux
    f1 = element vector due to flux boundary condition
    index1 = index for nodal dofs with flux
"""

# ---------------------------------
# input data for control parameters
# ---------------------------------
nel = 8            # number of elements
nnel = 4           # number of nodes per element
ndof = 1           # number of dofs per node
nnode = 15         # tota1 number of nodes m system
sdof = nnode*ndof  # total system dofs
deltt = 0.1        # time step size for transient analysis
stime = 0.0        # initial time
ftime = 1.0        # termination time
ntime = int(np.fix((ftime-stime) / deltt))  # number of time increment
a = 4266.7         # coefficient for the transient term
nf = 4             # number of element boundaries with flux
nnels = 2          # number of nodes per side of each element

# ---------------------------------------------
# input data for nodal coordinate values
# gcoord[i, j] where i->node No. and j-> x or y
# ---------------------------------------------
gcoord = np.array([
    [0.0,   0.0], [0.005,   0.0], [0.010,   0.0], [0.015,   0.0], [0.020,   0.0],
    [0.0, 0.005], [0.005, 0.005], [0.010, 0.005], [0.015, 0.005], [0.020, 0.005],
    [0.0,  0.01], [0.005,  0.01], [0.010,  0.01], [0.015,  0.01], [0.020,  0.01]])

# ---------------------------------------------------------
# input data for nodal connectivity for each element
# nodes[i, j] where i-> element No. and j-> connected nodes
# ---------------------------------------------------------
nodes = np.array([
    [0, 1,  6,  5], [1, 2,  7,  6], [2, 3,  8,  7], [3, 4,  9,  8],
    [5, 6, 11, 10], [6, 7, 12, 11], [7, 8, 13, 12], [8, 9, 14, 13]])

# ----------------------------------
# input data for boundary conditions
# ----------------------------------
# 1st, 5th, 6th, 10th, 11th and 15th nodes are constrained
bcdof = np.array([0, 4, 5, 9, 10, 14])  
# both described values are 100
bcval = np.array([300, 300, 300, 300, 300, 300])

# -------------------------------------------------------
# input for flux boundary conditions
# nflx[i, j] where i-> element no. and j-> two side nodes
# -------------------------------------------------------
# nodes on 1st, 2nd, 3rd and 4th element sides with flux
nflx = np.array([[10, 11], [11, 12], [12, 13], [13, 14]])
b = 100
c = 50

# --------------------------------------
# initialization of matrices and vectors
# --------------------------------------
ff = np.zeros(sdof)           # initialization of system vector
fn = np.zeros(sdof)           # initialization of effective system vector
fsol = np.zeros(sdof)         # solution vector
# time history solution of a selected node
sol = np.zeros(ntime+1, dtype='float64')       
kk = np.zeros([sdof, sdof])   # initialization of system matrix
mm = np.zeros([sdof, sdof])   # initialization of system matrix
index = np.zeros(nnel*ndof)   # initialization of index vector
f1 = np.zeros(nnels*ndof)               # element flux vector
k1 = np.zeros([nnels*ndof, nnels*ndof])      # flux matrix
index1 = np.zeros(nnels*ndof)           # flux index vector

# --------------------------------------------------------------
# computation of element matrices and vectors and their assembly
# --------------------------------------------------------------
nd = np.zeros(4, dtype='int64')
x = np.zeros(4)
y = np.zeros(4)
for iel in range(nel):
    for i in range(nnel):
        nd[i] = nodes[iel, i]   # connected nodes in (iel)-th element
        x[i] = gcoord[nd[i], 0] # x-coordinates of nodes in (iel)-th element
        y[i] = gcoord[nd[i], 1] # y-coordinates of nodes in (iel)-th element
    # extract system dofs for the element
    index = feeldof.feeldof(nd, nnel, ndof) 
    k = felp2dr4.felp2dr4(x, y)  # compute element matrix
    m = a * felpt2r4.felpt2r4(x, y)  # compute element matrix
    kk = feasmbl1.feasmbl1(kk, k, index) # assemble element matrices
    mm = feasmbl1.feasmbl1(mm, m, index) # assemble element matrices

# -----------------------------------------------------
# additional computation due to flux boundary condition
# -----------------------------------------------------
nds = np.zeros(2, dtype='int64')
for ifx in range(nf):
    nds[0] = nodes[ifx, 0]  # node with flux BC for (ifx)-th element
    nds[1] = nodes[ifx, 1]  # node with flux BC for (ifx)-th element
    x1 = gcoord[nds[0], 0]  # node coordinates
    y1 = gcoord[nds[0], 1]  # node coordinates
    x2 = gcoord[nds[1], 0]  # node coordinates
    y2 = gcoord[nds[1], 1]  # node coordinates
    eleng = np.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1))  # element side length
    index1 = feeldof.feeldof(nds, nnels, ndof) # find related system dofs
    k1 = b * feflxl2.feflxl2(eleng)     # compute e ement matrix due to flux
    f1 = b * c * fef1l.fef1l(eleng)  # compute element vector due to flux
    [kk,ff] = feasmbl2.feasmbl2(kk, ff, k1, f1, index1)  # assembly

# -------------------------
# loop for time integration
# -------------------------
# for ini in range(sdof):
#     fsol(in) = 0.0        # initial condition
fsol = np.ones(sdof) * 300
sol[0] = fsol[7]   # store time history solution for node 7
kk = mm + deltt*kk

for it in range(ntime):     # start loop for time integration
    fn = deltt*ff + np.dot(mm, fsol) # compute effective column vector
    kk, fn = feaplyc2.feaplyc2(kk, fn, bcdof, bcval) # apply boundary condition
    fsol = np.linalg.solve(kk, fn)  # solve the matrix equation 
    sol[it+1] = fsol[7]   # store time history solution for node 7
    
# ----------------------------
# plot the solution at nodes 7
# ----------------------------
time = np.arange(0, (ntime*deltt)+deltt, deltt)
fig, axis = plt.subplots()
axis.plot(time, sol, '*')
axis.set_xlabel('Time')
axis.set_ylabel('Solution at nodes 7')