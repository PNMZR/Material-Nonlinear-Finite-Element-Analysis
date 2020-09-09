# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 10:53:51 2020

@author: lenovo
"""
import numpy as np
import shape_func as sf
def Bmatrix4(pt,elemType,iel,element,node):
    """
    Gives the strain displacement matrix (B matrix of size 4x8) of each element
    global node element
    """

    sctr = element[iel-1, :]
    nn   = len(sctr)
    [N,dNdxi] = sf.shape_func(elemType,pt)  # element shape functions
    J0 = node[sctr-1,:].T@dNdxi             # element Jacobian matrix
    invJ0 = np.linalg.inv(J0)
    dNdx  = dNdxi@invJ0                  # derivatives of N w.r.t XY
    # Gpt = N.T@node(sctr-1,1)              # GP in global coord, used

    Bfem4 = np.zeros((4,2*nn))
    Bfem4[0,0:2*nn-1:2]  = dNdx[:,0] 
    Bfem4[1,1:2*nn:2]  = dNdx[:,1] 
    Bfem4[2,0:2*nn-1:2]  = dNdx[:,1] 
    Bfem4[2,1:2*nn:2]  = dNdx[:,0] 
    return Bfem4


