# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 21:22:29 2020

@author: lenovo
"""

import numpy as np
import gauss_rule as gr
import shape_func as sf

def selfwt_matrix(elemType, normal_order, gamma,node, element):
    """
    Generates the force matrix due to self weight
    """
    numelem = element.shape[0]
    numnode = node.shape[0]
    total_unknown = numnode*2
    selfwt = np.zeros(total_unknown)

    for iel in range(1, numelem+1):
        sctr1 = element[iel-1,:]      # element connectivity    
        swpt = sctr1*2                # element degree of freedom
   
        W, Q = gr.gauss_rule(iel, elemType, normal_order, element)                      
    
        for q in range(1, len(W)+1):
            pt = Q[q-1,:]
            wt = W[q-1]                             # quadrature point
            N, dNdxi = sf.shape_func(elemType,pt)
            J0 = node[sctr1-1,:].T@dNdxi          # element Jacobian matrix
            selfwt[swpt-1] = selfwt[swpt-1] + N*(-1*gamma)*np.linalg.det(J0)*wt # 都是在等参元里完成的
    return selfwt
