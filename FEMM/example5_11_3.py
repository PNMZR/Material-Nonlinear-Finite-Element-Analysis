# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import feeldof
import felp2dt3
import felpt2t3
import feasmbl1
import feaplyc2

"""
EX5.11.3.py
To solve the transient two-dimensional Laplace's equation using linear
triangular elements and backward difference method, given as:
    
           u,t = u,xx + u,yy  0 < x < 5, 0 < y < 2

Boundary conditions:
    
           u(0,y,t) = 100, u(5,y,t) = 100
           u,y(x,0,t) = 0, u,y(x,2,t) = 0
           
Initial condition:
           u(x,y,0) = 0 over the domain

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
"""

# ---------------------------------
# input data for control parameters
# ---------------------------------
nel = 16           # number of elements
nnel = 3           # number of nodes per element
ndof = 1           # number of dofs per node
nnode = 15         # tota1 number of nodes m system
sdof = nnode*ndof  # total system dofs
deltt = 0.4        # time step size for transient analysis
stime = 0.0        # initial time
ftime = 10         # termination time
ntime = int(np.fix((ftime-stime) / deltt))  # number of time increment

# ---------------------------------------------
# input data for nodal coordinate values
# gcoord[i, j] where i->node No. and j-> x or y
# ---------------------------------------------
gcoord = np.array([
    [0.0, 0.0], [1.25, 0.0], [2.5, 0.0], [3.75, 0.0], [5.0, 0.0],
    [0.0, 1.0], [1.25, 1.0], [2.5, 1.0], [3.75, 1.0], [5.0, 1.0],
    [0.0, 2.0], [1.25, 2.0], [2.5, 2.0], [3.75, 2.0], [5.0, 2.0]])

# ---------------------------------------------------------
# input data for nodal connectivity for each element
# nodes(i, j) where i-> element No. and j-> connected nodes
# ---------------------------------------------------------
nodes = np.array([
    [ 0,  1,  6], [ 1,  2,  7], [ 2,  3,  8], [ 3,  4,  9],
    [ 0,  6,  5], [ 1,  7,  6], [ 2,  8,  7], [ 3,  9,  8],
    [ 5,  6, 11], [ 6,  7, 12], [ 7,  8, 13], [ 8,  9, 14],
    [ 5, 11, 10], [ 6, 12, 11], [ 7, 13, 12], [ 8, 14, 13]])

# ----------------------------------
# input data for boundary conditions
# ----------------------------------
# 1st, 5th, 6th, 10th, 11th and 15th nodes are constrained
bcdof = np.array([0, 4, 5, 9, 10, 14])  
# both described values are 100
bcval = np.array([100, 100, 100, 100, 100, 100])

# --------------------------------------
# initialization of matrices and vectors
# --------------------------------------
ff = np.zeros(sdof)           # initialization of system vector
fn = np.zeros(sdof)           # initialization of effective system vector
fsol = np.zeros(sdof)         # solution vector
# vector containing time history solution
sol = np.zeros([2, ntime+1], dtype='float64')       
kk = np.zeros([sdof, sdof])   # initialization of system matrix
mm = np.zeros([sdof, sdof])   # initialization of system matrix
index = np.zeros(nnel*ndof)   # initialization of index vector


# --------------------------------------------------------------
# computation of element matrices and vectors and their assembly
# --------------------------------------------------------------
nd = np.zeros(3, dtype='int64')
x = np.zeros(3)
y = np.zeros(3)
for iel in range(nel):
    for i in range(nnel):
        nd[i] = nodes[iel, i]   # connected nodes in (iel)-th element
        x[i] = gcoord[nd[i], 0] # x-coordinates of nodes in (iel)-th element
        y[i] = gcoord[nd[i], 1] # y-coordinates of nodes in (iel)-th element
    # extract system dofs for the element
    index = feeldof.feeldof(nd, nnel, ndof) 
    k = felp2dt3.felp2dt3(x, y)  # compute element matrix
    m = felpt2t3.felpt2t3(x, y)  # compute element matrix
    kk = feasmbl1.feasmbl1(kk, k, index) # assemble element matrices
    mm = feasmbl1.feasmbl1(mm, m, index) # assemble element matrices
    
# -------------------------
# loop for time integration
# -------------------------
fsol = np.zeros(sdof)
sol[0, 0] = fsol[7]   # store time history solution for node 7
sol[1, 0] = fsol[8]   # store time history solution for node 8
kk = mm + deltt*kk

for it in range(ntime):     # start loop for time integration
    fn = deltt*ff + np.dot(mm, fsol) # compute effective column vector
    kk, fn = feaplyc2.feaplyc2(kk, fn, bcdof, bcval) # apply boundary condition
    fsol = np.linalg.solve(kk, fn)  # solve the matrix equation 
    sol[0, it+1] = fsol[7]   # store time history solution for node 7
    sol[1, it+1] = fsol[8]   # store time history solution for node 8

# ----------------------------------
# plot the solution at nodes 7 and 8
# ----------------------------------
time = np.arange(0, (ntime*deltt)+deltt, deltt)
fig, axis = plt.subplots()
axis.plot(time, sol[0, :], '*', label='Node 7')
axis.plot(time, sol[1, :], '-', label='Node 8')
axis.set_xlabel('Time')
axis.set_ylabel('Solution at nodes 7 and 8')
axis.legend()