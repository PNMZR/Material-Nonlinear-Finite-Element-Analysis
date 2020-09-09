# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 13:46:59 2020

@author: lenovo
"""
import gauss_pt_wt as gpw

def gauss_rule(iel, elemType, normal_order, element):
    """
    Provides the weight vector and gauss point coordinate matrix of
    the element based on the integration order selected.
    """

    #sctr = element[iel, :]                # element connectivity

    if ((elemType == 'Q4') and (normal_order <8)):
        W, Q = gpw.gauss_pt_wt(normal_order,'GAUSS',2)
    return W, Q

# function [W,Q] = gauss_rule(iel,elemType,normal_order)

# % Provides the weight vector and gauss point coordinate matrix of
# % the element based on the integration order selected.

# global node element

# sctr = element(iel,:);                % element connectivity


# if ((elemType == 'Q4') & (normal_order <8))
#     [W,Q] = gauss_pt_wt(normal_order,'GAUSS',2);
# elseif ((elemType == 'Q8') & (normal_order < 8))
#     [W,Q] = gauss_pt_wt(normal_order,'GAUSS',2);
# elseif ((elemType == 'Q9') & (normal_order < 8))
#     [W,Q] = gauss_pt_wt(normal_order,'GAUSS',2);
# elseif elemType == 'T3'
#     [W,Q] = gauss_pt_wt(normal_order,'TRIANGULAR',2);
# elseif elemType == 'T6'
#     [W,Q] = gauss_pt_wt(normal_order,'TRIANGULAR',2);
# end


# end  % end of function