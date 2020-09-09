# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 23:09:14 2020

@author: lenovo
"""
import numpy as np
import elementdof as eledof
import gauss_rule as gr
import shape_func as sf
import Bmatrix4 as bm4

def stress_calculation(node,element,elemType,U,normal_order,C):
    """
    calculates the element strains and stresses at the nodes
    in x, y and xy directions. 
    """
    numelem = element.shape[0]
    nonelm = element.shape[1]
    stress = np.zeros((4, nonelm, numelem))
    strain = np.zeros((4, nonelm, numelem))
    if elemType == 'Q4':
        stresspoints = np.array([[-1, -1], [ 1, -1],
                                  [ 1,  1], [-1,  1]])

    for iel in range(1, numelem+1):
        sctr = element[iel-1, :]   # element connectivity
        nn   = len(sctr)        # number of nodes per element
    
        eldof = eledof.elementdof(elemType, sctr)
        W, Q = gr.gauss_rule(iel, elemType, normal_order, element)   
        
        for kk in range(1, nn+1): 
            #pt = Q[kk-1,:]              # quadrature point 
            pt = stresspoints[kk-1, :]              # quadrature point      
            N, dNdxi = sf.shape_func(elemType, pt)  # element shape functions
            #J0 = node[sctr-1,:].T@dNdxi            # element Jacobian matrix       
            Bfem4 = bm4.Bmatrix4(pt, elemType, iel, element, node)           
            strain[:, kk-1, iel-1] = Bfem4@U[eldof-1]
            strain[3, :, :] = 0 
            stress[:, kk-1, iel-1] = C@strain[:, kk-1, iel-1]    
    return stress, strain

# function [stress,strain] =stress_calculation(node,element,elemType,U,normal_order,C)

# % calculates the element strains and stresses at the nodes
# % in x, y and xy directions. 

# numelem=size(element,1);nonelm=size(element,2);stress=zeros(4,nonelm,numelem);strain=zeros(4,nonelm,numelem);
# switch elemType
#     case 'Q9'
#       stresspoints=[-1 -1;1 -1;1 1;-1 1;0 -1;1 0;0 1;-1 0;0 0];
#     case 'Q8'
#       stresspoints=[-1 -1;1 -1;1 1;-1 1;0 -1;1 0;0 1;-1 0];
#     case 'Q4'
#       stresspoints=[-1 -1;1 -1;1 1;-1 1];
#     case 'T3'
#       stresspoints=[0 0;1 0;0 1];  
#     otherwise
#       stresspoints=[0 0;1 0;0 1;0.5 0;0.5 0.5;0 0.5];
# end
# for iel = 1 : numelem 
#     sctr = element(iel,:); % element connectivity
#     nn   = length(sctr);   % number of nodes per element
    
#    eldof =elementdof(elemType,sctr);  
#    [W,Q] = gauss_rule(iel,elemType,normal_order);   
        
#     for kk = 1:nn
#         pt = Q(kk,:);              % quadrature point      
#         [N,dNdxi] = shape_func(elemType,pt);  % element shape functions
#         J0 = node(sctr,:)'*dNdxi;             % element Jacobian matrix       
#         Bfem4 =Bmatrix4(pt,elemType,iel);           
#         strain(:,kk,iel)=Bfem4*U(eldof);
#         strain(4,:,:)=0; 
#         stress(:,kk,iel)=C*strain(:,kk,iel);  
      
#     end                  % end of looping on GPs
# end                      % end of looping on elements
# end   % end of function

