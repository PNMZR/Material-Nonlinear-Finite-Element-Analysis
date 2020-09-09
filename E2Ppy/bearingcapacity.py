# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 14:51:17 2020

@author: lenovo
"""

#  This script file is a FE code for a Rectangual element with 
#  4,8 or 9 nodes or for a triangular element with 3 or 6 nodes 
#  to do a bearing capacity problem for a Mohr Coulomb soil
# 节点编号和单元编号以及自由度编号均是从1开始，但是索引或者存储要减1

import time
import numpy as np
import matplotlib.pyplot as plt
import mesh_region as mesh
import supportcond as sup
import stiffness_matrix as sm
import selfwt_matrix as sem
import force_matrix as fm
import displacements as dis
import stress_calculation as sc
import invariants as invar
import elementdof as eledof
import gauss_rule as gr
import shape_func as sf
import Bmatrix4 as bm4
import invariants2 as invar2
import formm as formm
import formdg as formdg
import plot_m as plm
import plot_defo as pld
import plot_sig as pls

tic = time.perf_counter() 
elemType = 'Q4'
#--------------------------------------------------------------------
#                    INPUT PARAMETERS
#--------------------------------------------------------------------
E = 10e3          # 杨氏模量 KPa
nu = 0.3          # 泊松比
phi = 5           # 摩擦角 度
tsi = 5           # 剪胀角 度
c = 1             # 粘聚力 KPa
numx = 30          # x方向上的单元数量  节点数必须是11, 21, 31 ...
numy = 40          # y方向上的单元数量
D = 6             # 深度
L = 5             # 长度
gamma = 0         # 单位体积重量 KN/m3
df = -0.1         # 荷载增量
load_edge1 = 0    # 载荷起始点x坐标
load_edge2 = 1    # 载荷终止点x坐标
elemType = 'Q4'   # 单元类型
normal_order = 2  # 积分阶次
nsteps = 80       # 步数
maxit = 20        # 最大迭代数
tol = 1e-2        # 容许误差 
#--------------------------------------------------------------------

C = np.array([[1-nu,   nu,      0,   nu], 
              [  nu, 1-nu,      0,   nu], 
              [   0,    0, 0.5-nu,    0], 
              [   nu,  nu,      0, 1-nu]])*E/(1+nu)/(1-2*nu)
dt = 4*(1+nu)*(1-2*nu)/(E*(1-2*nu+np.sin(phi/5*np.pi)**2))
pt1 = np.array([0, -D/2])
pt2 = np.array([L, -D/2])
pt3 = np.array([L,  D/2])
pt4 = np.array([0,  D/2])

toc = time.perf_counter()
print(toc - tic,'seconds','    Mesh Generation....')

if elemType == 'Q4' or elemType == 'T3':
    node, element = mesh.mesh_region(pt1, pt2, pt3, pt4, numx, numy, elemType)

topEdge, dispNodes, dispNodes1, leftNodes1, topNodes = sup.supportcond(elemType, numx, numy) 

# Plot the FEM mesh 
plm.plot_m(elemType, dispNodes, dispNodes1, element, node, 'k-', deformed_mesh=False)
 
ko = 1 - np.sin(phi/180*np.pi) # k0 according to Jacky's equation
sigmatox = 0  # set horizontal loading zero 
numnode = node.shape[0]
numelem = element.shape[0]
nonelm = element.shape[1]
total_unknown = 2*numnode  
udofs  = np.hstack(((dispNodes*2)-1, (dispNodes1*2)-1)) # prescribed disp.in x-dir 位移已知自由度编号
vdofs  = dispNodes*2                             # prescribed disp. in y-dir
dofs = np.union1d(udofs[:],vdofs[:])                  # overall prescribed disp.位移已知自由度编号
unknowndof = np.setdiff1d(np.array(range(1,total_unknown+1)), dofs) # 自由度编号也是从1开始

toc1 = time.perf_counter()
print(toc1-toc,'seconds','     Stiffness Matrix Computation....')
K = sm.stiffness_matrix(node, element, elemType, normal_order, C)

toc2 = time.perf_counter()
print(toc2-toc1,'seconds','    Initial Strss with K0 Procedure....')
selfwt = sem.selfwt_matrix(elemType,normal_order,gamma,node,element)
sigmatoy1 = 0                        
          
if elemType == 'Q4':
    f, sctry = fm.force_matrix(node, topEdge, sigmatoy1, sigmatox, 
                                load_edge1, load_edge2)

U, u_x, u_y = dis.displacements(dispNodes, dispNodes1, numnode, K, f, selfwt)
stress, strain = sc.stress_calculation(node, element, elemType,
                                        U, normal_order, C)
stress[0, :, :] = ko*stress[1, :, :]  # 不理解
stress[3, :, :] = ko*stress[1, :, :]  # 不理解

if elemType == 'Q4':
    ule = numelem - numx + 1 # the element at the left top corner,索引时记得减去1
    stp = 1 # 第1个积分点， 索引时记得减去1
strainP = np.zeros((4, nonelm, numelem))     # set parameters to zero
stress_tr = np.zeros((4, nonelm, numelem))
ui = np.zeros(total_unknown)
load = 0
force = np.zeros(total_unknown)
f_old = np.zeros(total_unknown)
r = np.zeros(total_unknown)
b = np.zeros(total_unknown)
du = np.zeros(total_unknown)
dgds = np.zeros((4, nonelm, numelem))
stressule = stress[:, stp-1, ule-1]
p, q, theta = invar.invariants(stressule)

# prepare space for plotting data
pvq = np.zeros((2,nsteps+1))
pvq[0, 0] = p
pvq[1, 0] = q
evsyy = np.zeros((2,nsteps+1))
evsyy[0, 0] = 0
evsyy[1, 0] = stress[1,stp,ule] 
epsvq = np.zeros((2,nsteps+1))
epsvq[0, 0] = 0
epsvq[1, 0] = q
fvu = np.zeros((2,nsteps+1))
fvu[0, 0] = 0
fvu[1, 0] = 0                              
sigmatoy = 0

# start load stepping
for steps in range(1, nsteps+1):
    stepno = steps
    err = 1 # 初始误差设置为1，是不是得先计算
    nit = 0
    sigmatoy = df
    if elemType == 'Q4':
        f, sctry = fm.force_matrix(node, topEdge, sigmatoy, sigmatox,
                                    load_edge1, load_edge2)
    r = np.zeros(total_unknown)  # 残差        
    DEPS_PLA = np.zeros((4, nonelm, numelem))  # 塑性应变增量
    Du = np.zeros(total_unknown)
    du_old = np.zeros(total_unknown)
    Deps = np.zeros((4, nonelm, numelem))     # 应变增量
    dsig_pla = np.zeros((4, nonelm, numelem)) 
    Dsig = np.zeros((4, nonelm, numelem))     # 应力增量
    deps_pla = np.zeros((4, nonelm, numelem)) # 塑性应变率
   
    # start iteration loop
    while (err > tol) and (nit < maxit):
        nit = nit + 1
        du[unknowndof-1] = np.linalg.inv(K[np.ix_(unknowndof-1,unknowndof-1)])@f[unknowndof-1]         
        Du = Du + du
        for iel in range(1, numelem+1):           # start looping on all elements
            sctr = element[iel-1, :]              # element connectivity
            nn   = len(sctr)                      # number of nodes per element
            eldof = eledof.elementdof(elemType, sctr) 
            W, Q = gr.gauss_rule(iel, elemType, normal_order, element)     
            for kk in range(1, nn+1):  # start looping on all gauss points
                pt = Q[kk-1, :]                            
                N, dNdxi = sf.shape_func(elemType, pt) 
                J0 = node[sctr-1,:].T@dNdxi                   
                Bfem4 = bm4.Bmatrix4(pt, elemType, iel, element, node)
                Deps[:,kk-1,iel-1] = Bfem4@du[eldof-1]  # 应变增量
                Deps[-1,-1,:] = np.zeros(Deps.shape[2])
                Dsig[:,kk-1,iel-1] = C@(Deps[:,kk-1,iel-1] - DEPS_PLA[:,kk-1,iel-1]) # 为什么减去
                stress_tr[:,kk-1,iel-1] = stress[:,kk-1,iel-1] + Dsig[:,kk-1,iel-1]              
                p, q, theta = invar2.invariants2(kk,iel,stress_tr)
                F = (p*np.sin(phi/180*np.pi) 
                      + q*((np.cos(theta/180*np.pi)/np.sqrt(3))
                      - (np.sin(theta/180*np.pi)*np.sin(phi/180*np.pi)/3))
                      - c*np.cos(phi/180*np.pi))
                if F < 0:
                    stress[:,kk-1,iel-1] = stress_tr[:,kk-1,iel-1]
                    err = 0 # 为啥误差就设置为0了
                else:     
                    m1,m2,m3 = formm.formm(kk,iel,stress_tr)
                    dg1, dg2, dg3 = formdg.formdg(tsi,q,theta)
                    dgds[:, kk-1, iel-1] = (dg1*m1 + dg2*m2 + dg3*m3)@stress_tr[:,kk-1,iel-1]
                    deps_pla[:,kk-1,iel-1] = dt*F*dgds[:,kk-1,iel-1]         
                    DEPS_PLA[:,kk-1,iel-1] = DEPS_PLA[:,kk-1,iel-1] + deps_pla[:,kk-1,iel-1]  
                    if nit == 1:
                        err = 1
                    else:
                        err = np.max(np.abs(du_old[2:-1:2]-du[2:-1:2]))  
                r[eldof-1] = r[eldof-1]+Bfem4.T@(C@deps_pla[:,kk-1,iel-1])*W[kk-1]*np.linalg.det(J0)          
            f[sctry-1] = f[sctry-1] + r[sctry-1]
        du_old = du
        f_old = f
    
    stress[:,kk-1,iel-1] = stress_tr[:,kk-1,iel-1]      
    s = np.transpose(stress,[1, 0, 2])
    sigma = np.transpose(stress,[2, 1, 0])  
    ui = ui + Du
    u_x = ui[0:2*numnode-1:2]
    u_y = ui[1:2*numnode:2]
    strain = strain + Deps    
    strainP = strainP + DEPS_PLA    
    xxx = strain[1, stp-1, ule-1]
    yyy = stress[1, stp-1, ule-1]               
    eps_vol = strainP[0, stp-1, ule-1] + strainP[1, stp-1, ule-1] + strainP[3, stp-1, ule-1]
    stressule = stress[:, stp-1, ule-1] 
    p, q, theta = invar.invariants(stressule)
    pvq[0, steps+1-1] = p          # 平均应力
    pvq[1, steps+1-1] = q          # 偏应力
    epsvq[0, steps+1-1] = eps_vol  # 体积应变
    epsvq[1, steps+1-1] = q        # 偏应力 
    evsyy[0, steps+1-1] = xxx      # epsilon_y
    evsyy[1, steps+1-1] = yyy      # sigma_y
    load = load + df 
    fvu[0, steps+1-1] = u_y[-1]
    fvu[1, steps+1-1] = load

fig1, ax1 = plt.subplots()
ax1.plot(abs(pvq[0, :]), pvq[1, :], '--b*', linewidth=2)

ax1.set_xlabel('P', fontsize=16);
ax1.set_ylabel('q', fontsize=16);

fig2, ax2 = plt.subplots()
ax2.plot(abs(epsvq[0, :]), epsvq[1, :], '--r*', linewidth=2)

ax2.set_xlabel('epsvol', fontsize=16)
ax2.set_ylabel('q', fontsize=16)

fig3, ax3 = plt.subplots()
ax3.plot(abs(evsyy[0, :]), abs(evsyy[1, :]), '--r*', linewidth=2)

ax3.set_xlabel('epsyy', fontsize=16)
ax3.set_ylabel('sigyy', fontsize=16)  


fig4, ax4 = plt.subplots()
ax4.plot(abs(fvu[0, :]), abs(fvu[1, :]), '--r*', linewidth=2)

ax4.set_xlabel('u_y', fontsize=16)
ax4.set_ylabel('load', fontsize=16)  

# Plot numerical deformed configuration
dispnorm = L/np.max(np.sqrt(u_x**2+u_y**2))
fac = dispnorm*0.05   # magnification factor
plm.plot_m(elemType, dispNodes, dispNodes1, element,
                node+fac*np.vstack([u_x, u_y]).T, 'k-', deformed_mesh=True)

# plot stress and deformation intensity with a colormap
pld.plot_defo(pt1, pt2, pt3, pt4, numx, numy, fac, u_x, u_y)
pls.plot_sig(pt1, pt2, pt3, pt4, numx, numy, fac, u_x, u_y, sigma, element)