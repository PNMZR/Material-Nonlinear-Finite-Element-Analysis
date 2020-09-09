""""
EX3.5.3.py
To solve the ordinary differential equation given as
x^2*u" + 2 * x * u'+ 4 * u = x^2, 10 < x < 20
u(1O) = 0 and u(20) = 100
using 10 linear elements

Variable descriptions
k = element matrix
f = element vector
kk = system matnx
ff = system vector
index = a vector containing system dofs associated with each element
bcdof = a vector containing dofs associated with boundary conditions
bcval = a vector containing boundary condition values associated with
        the dofs in bcdof
"""

import numpy as np
import feaplyc2
import feeldof1
import feodex2l
import fefx2l
import feasmbl2
import matplotlib.pyplot as plt

# ---------------------------------
# input data for control parameters
# ---------------------------------
nel10 = 10              # number of elements
nnel = 2                # number of nodes per element
ndof = 1                # number of dofs per node
nnode10 = 11            # total number of nodes in system
sdof = nnode10 * ndof   # total system dofs

# --------------------------------------
# input data for nodal coordinate values
# --------------------------------------

gcoord_10 = np.array([10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])

# --------------------------------------------------
# input data for nodal connectivity for each element
# --------------------------------------------------
nodes_10 = np.array([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8],
                     [8, 9], [9, 10], [10, 11]])

# ----------------------------------
# input data for boundary conditions
# ----------------------------------
bcdof = np.array([1, 11])  # first and 11th nodes are constrained
bcval = np.array([0, 100])  # whose described value is 0

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
for iel in range(nel10):                # loop for the total number of elements
    nd[0] = nodes_10[iel, 0]                  # extract nodes for (iel)-th element
    nd[1] = nodes_10[iel, 1]  
    xl = gcoord_10[nd[0]-1]                     # extract nodal coord values
    xr = gcoord_10[nd[1]-1]
    eleng = xr - xl                     #  eleng - element length
    index = feeldof1.feeldof1(nd, nnel, ndof)   # extract system dofs associated

    k = feodex2l.feodex2l(xl, xr, eleng)      # compute element matrix
    f = fefx2l.fefx2l(xl, xr, eleng)    # compute element vector
    # assemble element matrices and vectors
    kk, ff = feasmbl2.feasmbl2(kk, ff, k, f, index)
   
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
    esol_i = 0.00102 * x ** 4 - 0.16667 * x ** 2 + 64.5187 / x
    eso.append(esol_i)
eso = np.array(eso)

# ----------------------------------
# print both exact and fem solutions
# ----------------------------------
num = range(1, 1, sdof+1)

fig, axis = plt.subplots()
axis.plot(gcoord_10, eso, 'r--', label="Analytical solution", linewidth=4) 
axis.plot(gcoord_10, fsol, 'k-', label="FEM solution")
plt.rcParams['font.sans-serif'] = ['SimSun']
 # 解决保存图像是负号'-'显示为方块的问题
plt.rcParams['axes.unicode_minus'] = False 
axis.set_xlabel("坐标 $x$")
axis.set_ylabel("$u$")
axis.legend()