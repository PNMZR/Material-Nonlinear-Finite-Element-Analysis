import numpy as np

def feeldof(nd, nnel, ndof):
    """-------------------------------------------------------------
    Purpose:
    Compute system dofs associated with each element

    Synopsis:
    index = feeldof(nd, nnel, ndof)

    Variable Description:
    index - system dof vector associated with element iel
    nd - element node numbers whose system dofs are to be determined
    nnel - number of nodes per element
    ndof - number of dofs per node
    -------------------------------------------------------------"""
    index = []
    for i in range(nnel):
        start = (nd[i]-1)*ndof
        for j in range(1, ndof+1):
            index_j = start + j
            index.append(index_j)
    return np.array(index)