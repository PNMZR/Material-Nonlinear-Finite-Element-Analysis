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
EX5.11.5.py
To solve the transient two-dimensional Laplace's equation using bilinear
rectangular elements and the Crank-Nicolson technique, given as:
    
           0.04*u,t = u,xx + u,yy  0 < x < 5, 0 < y < 10
           
Boundary conditions:
    
           u(x,0) = 0, u(x,10) = lOOsin(pi*x/10),
           u,y(x,0,t) = 0, u,y(x,0.01,t) = 20(u-50)
           
Initial condition:
    
           u(x,y,0) = 100 over the domain

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
nnel = 4           # number of nodes per element
ndof = 1           # number of dofs per node
nnode = 25         # tota1 number of nodes m system
sdof = nnode*ndof  # total system dofs
deltt = 0.04       # time step size for transient analysis
stime = 0.0        # initial time
ftime = 2          # termination time
ntime = int(np.fix((ftime-stime) / deltt))  # number of time increment
a = 0.04         # coefficient for the transient term

# ---------------------------------------------
# input data for nodal coordinate values
# gcoord[i, j] where i->node No. and j-> x or y
# ---------------------------------------------
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
kn = np.zeros([sdof, sdof])   # effective system matrix
index = np.zeros(nnel*ndof)   # initialization of index vector

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

# -------------------------
# loop for time integration
# -------------------------
fsol = np.ones(sdof) * 100.0
sol[0] = fsol[12]     # store time history solution for node 12
kn = 2*mm + deltt*kk  # compute effective system matrix

for it in range(ntime):     # start loop for time integration
    # compute effective column vector
    fn = deltt*ff + np.dot((2*mm - deltt*kk), fsol) 
    kn, fn = feaplyc2.feaplyc2(kn, fn, bcdof, bcval) # apply boundary condition
    fsol = np.linalg.solve(kn, fn)  # solve the matrix equation 
    sol[it+1] = fsol[12]   # store time history solution for node 12
    
# ----------------------------
# plot the solution at nodes 7
# ----------------------------
time = np.arange(0, (ntime*deltt)+deltt, deltt)
fig, axis = plt.subplots()
axis.plot(time, sol, '*')
axis.set_xlabel('Time')
axis.set_ylabel('Solution at nodes 7')