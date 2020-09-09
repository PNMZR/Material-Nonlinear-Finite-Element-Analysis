# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 15:06:01 2020

@author: lenovo
"""
import numpy as np
import shape_func as sf

def Bmatrix(pt, elemType, iel, element, node):
    """
    Gives the strain displacement matrix (B matrix of size 3x8)of each element
    global node element
    """

    sctr = element[iel-1,:]
    nn   = len(sctr)
    N, dNdxi = sf.shape_func(elemType, pt)   # element shape functions
    J0 = node[sctr-1,:].T@dNdxi                # element Jacobian matrix
    invJ0 = np.linalg.inv(J0)
    dNdx  = dNdxi@invJ0                      # derivatives of N w.r.t XY
    #Gpt = N.T@node[sctr,:]                   # GP in global coord, used

    Bfem = np.zeros((3,2*nn))
    Bfem[0, 0:-1:2] = dNdx[:,0].T
    Bfem[1, 1::2] = dNdx[:,1].T
    Bfem[2, 0:-1:2] = dNdx[:,1].T
    Bfem[2, 1::2] = dNdx[:,0].T

    return Bfem  # 应变位移矩阵   

# function Bfem =Bmatrix(pt,elemType,iel)
# % Gives the strain displacement matrix (B matrix of size 3x8)of each element
# global node element

# sctr = element(iel,:);
# nn   = length(sctr);
# [N,dNdxi] = shape_func(elemType,pt);  % element shape functions
# J0 = node(sctr,:)'*dNdxi;             % element Jacobian matrix
# invJ0 = inv(J0);
# dNdx  = dNdxi*invJ0;                  % derivatives of N w.r.t XY
# Gpt = N'*node(sctr,:);               % GP in global coord, used


#     Bfem = zeros(3,2*nn);
#     Bfem(1,1:2:2*nn)  = dNdx(:,1)' ;
#     Bfem(2,2:2:2*nn)  = dNdx(:,2)' ;
#     Bfem(3,1:2:2*nn)  = dNdx(:,2)' ;
#     Bfem(3,2:2:2*nn)  = dNdx(:,1)' ;

# end  % end of function                    


    
    