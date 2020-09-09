# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 09:22:44 2020

@author: lenovo
"""
import numpy as np
import elementdof as eledof
import gauss_rule as gr
import shape_func as sf
import Bmatrix as bm

def stiffness_matrix(node,element,elemType,normal_order,C):
    """
    Generates the element and global stiffness matrix
    """
    numnode = node.shape[0]
    numelem = element.shape[0]
    total_unknown = numnode*2            
    K = np.zeros((total_unknown,total_unknown))          

    for iel in range(1, numelem+1): 
        sctr = element[iel-1, :]                   # element connectivity   
        eldof = eledof.elementdof(elemType,sctr)   # element degree of freedom
    
        W, Q = gr.gauss_rule(iel,elemType,normal_order,element)  # 这块感觉ie和element是多余的           
      
        
        for kk in range(1, len(W)+1):          # Loop on Gauss points 积分节点个数也是从1开始，索引时减去1
            pt = Q[kk-1,:]                  # quadrature point
            N, dNdxi = sf.shape_func(elemType, pt)  
            J0 = node[sctr-1,:].T@dNdxi               # element Jacobian matrix
            Bfem = bm.Bmatrix(pt,elemType,iel,element,node)
       
            K[np.ix_(eldof-1,eldof-1)] += Bfem.T@C[0:3,0:3]@Bfem*W[kk-1]*np.linalg.det(J0)
    return K

