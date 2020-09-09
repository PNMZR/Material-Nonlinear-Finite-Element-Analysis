# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 23:06:38 2020

@author: lenovo
"""
import numpy as np

def displacements(dispNodes, dispNodes1, numnode, K, f, selfwt):
    """
    evaluates the unknown degree of freedom (displacements) at the nodes
    """
    total_unknown = 2*numnode
    udofs  = np.hstack(((dispNodes*2)-1,(dispNodes1*2)-1)) # prescribed disp.in x-dir
    vdofs  = dispNodes*2                             # prescribed disp. in y-dir
    dofs = np.union1d(udofs[:],vdofs[:])                  # overall prescribed disp.
    unknowndof = np.setdiff1d(np.array(range(1,total_unknown+1)), dofs)
    
    

    F = f[unknowndof-1] + selfwt[unknowndof-1]
    u = np.linalg.inv(K[np.ix_(unknowndof-1, unknowndof-1)])@F
    U = np.zeros(total_unknown)
    U[unknowndof-1] = u
    u_x = U[0:2*numnode-1:2]        # 0 2 4 8 ...2*(numnode-1)
    u_y = U[1:2*numnode:2]          # 1 3 5 7 ...2*numnode-1实际的索引

    return U, u_x, u_y

# function [U,u_x,u_y] =displacements(dispNodes,dispNodes1,numnode,K,f,selfwt)

# % evaluates the unknown degree of freedom (displacements) at the nodes

# total_unknown=2*numnode;

# udofs  = [(dispNodes.*2)-1;(dispNodes1.*2)-1]; %prescribed disp.in x-dir
# vdofs  = dispNodes.*2;                         %prescribed disp. in y-dir
# dofs=union(udofs(:),vdofs(:));                 %overall prescribed disp.
# unknowndof=setdiff((1:total_unknown)',dofs);

# F=f(unknowndof)+selfwt(unknowndof);
# u=K(unknowndof,unknowndof)\F;
# U=zeros(total_unknown,1);
# U(unknowndof)=u;
# u_x = U(1:2:2*numnode-1) ; % 1 3 5 7 ...
# u_y = U(2:2:2*numnode) ; % 2 4 6 8 ...

# end    % end of function

