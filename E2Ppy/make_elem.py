# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 18:04:06 2020

@author: lenovo
"""
import numpy as np

def make_elem(node_pattern, numx, numy, inc_u, inc_v):
    """
    creates a connectivity list of primary nodes in Q4 and T3 element
    """
    inc = 0
    e = 1
    element = np.zeros((numx*numy, len(node_pattern)), dtype=np.int64)

    for row in range(1, numy+1):
        for col in range(1, numx+1):
            element[e-1, :] = node_pattern + inc
            inc = inc + inc_u  # 下一次执行的inc
            e = e + 1
        inc = row*inc_v
    return element

# function element = make_elem(node_pattern, numx, numy, inc_u, inc_v)

# % creates a connectivity list of primary nodes in Q4 and T3 element

# if ( nargin < 5 )
#     disp('Not enough parameters specified for make_elem function')
# end

# inc = [zeros(1, size(node_pattern,2))];
# e = 1;
# element = zeros(numx*numy, size(node_pattern,2));

# for row = 1 : numy
#     for col = 1 : numx
#         element(e, :) = node_pattern + inc;
#         inc = inc + inc_u;
#         e = e + 1;
#     end
#     inc = row * inc_v;
# end
# end
