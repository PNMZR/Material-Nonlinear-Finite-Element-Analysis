# -*- coding: utf-8 -*-

def feasmbl2(kk, ff, k, f, index):
    """
    Purpose:
        Assembly of element matrices into the system matrix and assembly of
        element vectors into the system vector.
    
        Synopsis: kk, ff = feasmbl2(kk, ff ,k,f,index)

        Variable Description:
            kk - system matrix
            ff - system vector
            k - element matrix
            f - element vector
            index - d.o.f. vector associated with an element
    """

    edof = len(index)
    for i in range(edof):
        ii = index[i]
        ff[ii] = ff[ii] + f[i]
        for j in range(edof):
            jj = index[j]
            kk[ii, jj] = kk[ii, jj] + k[i, j]        
    return kk, ff