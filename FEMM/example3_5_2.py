#----------------------------------------------------------------------
# EX3.5.2.py
# To solve the ordinary differential equation given as
# a u" + b u'+ c u = 1, 0 < x < 1
# u(O) = 0 and u'(1) = 0
# using 5 or 10 linear elements
#
# Variable descriptions
# k = element matrix
# f = element vector
# kk = system matnx
# ff = system vector
# index = a vector containing system dofs associated with each element
# bcdof = a vector containing dofs associated with boundary conditions
# bcval = a vector containing boundary condition values associated with
# the dofs in bcdof
# ---------------------------------------------------------------------

import numpy as np
import feaplyc2
import feeldof1
import feode2l
import fef1l
import feasmbl2
import matplotlib.pyplot as plt

# ---------------------------------
# input data for control parameters
# ---------------------------------
nel5 = 5              # number of elements
nel10 = 10              # number of elements
nnel = 2             # number of nodes per element
ndof = 1             # number of dofs per node
nnode5 = 6
nnode10 = 11            # total number of nodes in system
sdof = nnode10 * ndof  # total system dofs

# --------------------------------------
# input data for nodal coordinate values
# --------------------------------------
gcoord_5 = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

gcoord_10 = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

# --------------------------------------------------
# input data for nodal connectivity for each element
# --------------------------------------------------
nodes_5 = np.array([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6]])
nodes_10 = np.array([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8],
                     [8, 9], [9, 10], [10, 11]])

# --------------------------------------
# input data for coefficients of the ODE
# --------------------------------------
acoef = 1     # coefficient 'a' of the diff eqn
bcoef = -3    # coefficient 'b' of the diff eqn
ccoef = 2     # coefficient 'c' of the diff eqn

# ----------------------------------
# input data for boundary conditions
# ----------------------------------
bcdof = np.array([1])  # first is constrained
bcval = np.array([0])  # whose described value is 0

# --------------------------------------
# initialization of matrices and vectors
# --------------------------------------
ff = np.zeros([sdof, 1])            # initialization of system force vector
kk = np.zeros([sdof, sdof])         # initialization of system matrix
index = np.zeros([nnel*ndof, 1])    # initialization of index vector  

# --------------------------------------------------------------
# computation of element matrices and vectors and their assembly
# --------------------------------------------------------------
nd = np.zeros(2, dtype='int64')
for iel in range(nel10):
    nd[0] = nodes_10[iel, 0]                  # extract nodes for (iel)-th element
    nd[1] = nodes_10[iel, 1]  
    xl = gcoord_10[nd[0]-1]                     # extract nodal coord values
    xr = gcoord_10[nd[1]-1]
    eleng = xr - xl                     #  eleng - element length
    index = feeldof1.feeldof1(nd, nnel, ndof)   # extract system dofs associated

    k = feode2l.feode2l(acoef, bcoef, ccoef, eleng)  # compute element matrix
    f = fef1l.fef1l(xl, xr, eleng)                        # compute element vector
    # assemble element matrices and vectors
    kk, ff = feasmbl2.feasmbl2(kk, ff, k, f, index)

#------------------------------------------------------
# apply the natural boundary condition at the last node
# -----------------------------------------------------
ff[nnode10-1] = ff[nnode10-1] - 1
   
# -------------------------
# apply boundary conditions
# -------------------------

kk, ff = feaplyc2.feaplyc2(kk, ff, bcdof , bcval)

# -------------------------
# solve the matrix equation
# -------------------------
fsol = np.linalg.solve(kk, ff)

# -------------------
# analytical solution
# -------------------
cl = (1 + 0.5 * np.exp(1)) / (2 * np.exp(2) - np.exp(1))
c2 = -(1 + np.exp(2)) / (2 * np.exp(2) - np.exp(1))
eso = []
for i in range(nnode10):
    x = gcoord_10[i]
    esol_i = cl * np.exp(2 * x) + c2 * np.exp(x) + 1 / 2
    eso.append(esol_i)
eso = np.array(eso)

# ----------------------------------
# print both exact and fem solutions
# ----------------------------------
num = range(1, 1, sdof+1)

fig, axis = plt.subplots()
axis.plot(gcoord_10, eso, 'r--', label="Analytical solution", linewidth=4) # element modified
axis.plot(gcoord_10, fsol, 'k-', label="FEM solution")        # element modified
plt.rcParams['font.sans-serif'] = ['SimSun']
 # 解决保存图像是负号'-'显示为方块的问题
plt.rcParams['axes.unicode_minus'] = False 
axis.set_xlabel("坐标 $x$")
axis.set_ylabel("$u$")
axis.legend()