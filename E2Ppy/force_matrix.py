# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 22:16:44 2020

@author: lenovo
"""
import numpy as np
import gauss_pt_wt as gpw
import shape_func as sf

def force_matrix(node,topEdge,sigmatoy,sigmatox,load_edge1,load_edge2):

    """
    Generates the force matrix due to externally applied loads
    for Q4 and T3 elements and the location of these loads
    """

    numnode = node.shape[0]
    total_unknown = numnode*2
    f = np.zeros(total_unknown)

    W, Q = gpw.gauss_pt_wt(1,'GAUSS',1);
 
    ii1 = np.intersect1d(np.where(node[:, 0] == load_edge1)[0]+1,np.unique(topEdge))
    jj1 = np.intersect1d(np.where(node[:, 0] == load_edge2)[0]+1,np.unique(topEdge))
    ii2 = np.where(topEdge[:, 0] == ii1)[0]+1
    jj2 = np.where(topEdge[:, 1] == jj1)[0]+1

    for e in range(int(ii2),int(jj2)+1):
        sctr = topEdge[e-1,:]
        sctry = sctr*2
        sctrx = sctr*2-1
        for q in range(len(W)):
            pt = Q[q,:]
            wt = W[q]
            N, dNdxi  = sf.shape_func('L2', pt)
            J0 = np.abs(node[sctr[1],0]-node[sctr[0], 0])/2
            f[sctry-1] = f[sctry-1] + N*sigmatoy*J0*wt
            f[sctrx-1] = f[sctrx-1] + N*sigmatox*J0*wt
    return f, sctry   # 不太好理解，换另一种方式p236 fish book 为什么要输出sctry

# function [f,sctry]=force_matrix(node,topEdge,sigmatoy,sigmatox,load_edge1,load_edge2)

# % Generates the force matrix due to externally applied loads
# % for Q4 and T3 elements and the location of these loads

# numnode = size(node,1);
# total_unknown = numnode*2;
# f = zeros(total_unknown,1);

# [W,Q]=gauss_pt_wt(1,'GAUSS',1);

# ii1=intersect(find(node(:,1)==load_edge1),unique(topEdge));
# jj1=intersect(find(node(:,1)==load_edge2),unique(topEdge));
# ii2=find(topEdge(:,1)==ii1);
# jj2=find(topEdge(:,2)==jj1);

# for e =ii2:jj2
#        sctr = topEdge(e,:);
#        sctry = sctr.*2 ;
#        sctrx = sctr.*2-1;
#     for q=1:size(W,1)
#         pt = Q(q,:);
#         wt = W(q);
#         N  = shape_func('L2',pt);
#         J0 = abs( node(sctr(2))-node(sctr(1)) )/2;
#         f(sctry) = f(sctry) + N*sigmatoy*det(J0)*wt;
#         f(sctrx) = f(sctrx) + N*sigmatox*det(J0)*wt;
#     end   % end of quadrature loop
# end       % end of element loop

# end  % end of function

