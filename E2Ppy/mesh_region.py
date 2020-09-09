# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 16:59:36 2020

@author: lenovo
"""
import numpy as np
import square_node_array as sna
import make_elem as me

def mesh_region(pt1, pt2, pt3, pt4, numx, numy, elemType):
    """
    Generates an array of nodal connectivity (coordinates of each node)
    and element connectivity (nodes of each element with counterclockwise
    ordering.)
    given the four corners points of the domain, number of elements
    in each direction (numx, numy) and the element type (Q4,Q8,Q9 and T3)
    """

    nnx = numx + 1
    nny = numy + 1

    if elemType == 'Q4':
        node = sna.square_node_array(pt1, pt2, pt3, pt4, nnx, nny)
        inc_u = 1
        inc_v = nnx
        node_pattern = np.array([1, 2, nnx+2, nnx+1], dtype=np.int64)
        element = me.make_elem(node_pattern, numx, numy, inc_u, inc_v)    
    elif elemType == 'T3':
        node = sna.square_node_array(pt1, pt2, pt3, pt4, nnx, nny)
        node_pattern1 = np.array([1, 2, nnx+1], dtype=np.int64)
        node_pattern2 = np.array([2, nnx+2, nnx+1], dtype=np.int64)
        inc_u = 1
        inc_v = nnx
        numberelem = 2*numx*numy
        element = np.zeros((numberelem, 3), dtype=np.int64)
        element1 = me.make_elem(node_pattern1, numx,numy, inc_u, inc_v)
        element2 = me.make_elem(node_pattern2, numx,numy, inc_u, inc_v)
        element[0:-1:2,:] = element1
        element[1::2,:] = element2
    return node, element


# function [node,element] = mesh_region(pt1, pt2, pt3, pt4,numx,numy,elemType)

# % Generates an array of nodal connectivity (coordinates of each node)
# % and element connectivity (nodes of each element with counterclockwise
# % ordering.)
# % given the four corners points of the domain,number of elements
# % in each direction (numx,numy)and the element type (Q4,Q8,Q9 and T3)

# global L D
# nnx = numx + 1;
# nny = numy + 1;

# switch elemType
    
#     case 'Q4'
#         node = square_node_array(pt1, pt2, pt3, pt4, nnx, nny);
#         inc_u = 1;
#         inc_v = nnx;
#         node_pattern = [1 2 nnx+2 nnx+1];
#         [element] = make_elem(node_pattern, numx, numy, inc_u, inc_v);
        
#     case 'Q9'
#         [node,element]=structured_q9_mesh(pt1,pt2,pt3,pt4,numx,numy);
        
#     case 'Q8'
#         [node,element]=structured_q8_mesh(pt1,pt2,pt3,pt4,numx,numy);
        
#     case 'T3'
        
#         node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny);
#         node_pattern1=[ 1 2 nnx+1 ];
#         node_pattern2=[ 2 nnx+2 nnx+1 ];
#         inc_u=1;
#         inc_v=nnx;
#         numberelem=2*numx*numy;
#         element=zeros(numberelem,3);
#         element1=make_elem(node_pattern1,numx,numy,inc_u,inc_v);
#         element2=make_elem(node_pattern2,numx,numy,inc_u,inc_v);
#         element(1:2:end,:)=element1;
#         element(2:2:end,:)=element2;
#     otherwise
#         error('only Q4,Q9,Q8,and T3 are supported by this function');
# end
# end