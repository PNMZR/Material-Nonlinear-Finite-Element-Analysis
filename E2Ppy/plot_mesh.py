# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 16:12:24 2020

@author: lenovo
"""
import numpy as np
import matplotlib.pyplot as plt

def plot_mesh(node, element, elem_type, se):
    """
    Plots the finite element mesh 
    """

    # fill node if needed
    # if (size(node, 2) < 3)
    #    for c = size(node,2) + 1:3
    #       node(:, c)=[zeros(size(node,1), 1)]

    fig, ax = plt.subplots()
    ax.set_title('Undeformed FE mesh')
    for e in range(1, element.shape[0]+1):
        if elem_type == 'Q9':                   # 9-node quad element
            ord = [1, 5, 2, 6, 3, 7, 4, 8, 1]
        elif elem_type == 'Q8':                 # 8-node quad element
            ord = [1, 5, 2, 6, 3, 7, 4, 8, 1]
        elif elem_type == 'T3':                 # 3-node triangle element
            ord = [1, 2, 3, 1]
        elif elem_type == 'T6':                 # 6-node triangle element
            ord = [1, 4, 2, 5, 3, 6, 1]
        elif elem_type == 'Q4':                  # 4-node quadrilateral element
            ord = [1, 2, 3, 4, 1]
        elif elem_type == 'L2':                  # 2-node line element
            ord = [1, 2]   
  
    
        xpt = np.zeros(len(ord))
        ypt = np.zeros(len(ord))
        for n in range(1, len(ord)+1):
            xpt[n-1] = node[element[e-1, ord[n-1]-1]-1, 0]
            ypt[n-1] = node[element[e-1, ord[n-1]-1]-1, 1]     
        ax.plot(xpt, ypt, se)
      

# function plot_mesh(node,connect,elem_type,se)

# % Plots the finite element mesh 
  
# if ( nargin < 4 )
#    se='w-';
# end

# holdState=ishold;
# hold on

# % fill node if needed
# if (size(node,2) < 3)
#    for c=size(node,2)+1:3
#       node(:,c)=[zeros(size(node,1),1)];
#    end
# end

# for e=1:size(connect,1)
  
#    if ( strcmp(elem_type,'Q9') )       % 9-node quad element
#       ord=[1,5,2,6,3,7,4,8,1];
#    elseif ( strcmp(elem_type,'Q8') )  % 8-node quad element
#       ord=[1,5,2,6,3,7,4,8,1];
#    elseif ( strcmp(elem_type,'T3') )  % 3-node triangle element
#       ord=[1,2,3,1];
#    elseif ( strcmp(elem_type,'T6') )  % 6-node triangle element
#       ord=[1,4,2,5,3,6,1];
#    elseif ( strcmp(elem_type,'Q4') )  % 4-node quadrilateral element
#       ord=[1,2,3,4,1];
#    elseif ( strcmp(elem_type,'L2') )  % 2-node line element
#       ord=[1,2];   
#    end
   
#    for n=1:size(ord,2)
#       xpt(n)=node(connect(e,ord(n)),1);
#       ypt(n)=node(connect(e,ord(n)),2);      
#       zpt(n)=node(connect(e,ord(n)),3);
#    end
#    plot3(xpt,ypt,zpt,se)
  
# end
#     axis equal
      
# if ( ~holdState )
#   hold off
# end
# end % end of function
