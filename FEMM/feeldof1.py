# -*- coding: utf-8 -*-

import numpy as np

def feeldof1(nd, nnel, ndof):
    """------------------------------------------------------------------
    Purpose:
    Compute system dofs associated with each element in one-dimensional
    problem

    Synopsis:
    index = feeldofl(iel, nnel, ndof)
    
    Variable Description:
    index - system dof vector associated with element iel
    iel - element number whose system dofs are to be determined
    nnel - number of nodes per element
    ndof - number of dofs per node
    ------------------------------------------------------------------"""
 
    index = []
    for i in range(nnel):
        start = (nd[i]-1)*ndof
        for j in range(ndof):
            index_j = start + j + 1
            index.append(index_j)
    return np.array(index)