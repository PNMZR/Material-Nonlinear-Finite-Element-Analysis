# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 10:51:22 2020

@author: lenovo
"""
import numpy as np

def supportcond(elemType, numx, numy):
    """
    Forms a vector of nodes which are restrained in x direction or in both
    x and y direction which are the supports of the domain. 
    It also forms a matrix of nodes for load application at the top of the 
    domain.
    """

    if elemType == 'Q4' or elemType == 'T3':
        nnx = numx + 1
        nny = numy + 1
        uln = nnx*(nny-1) + 1      # upper left node number
        urn = nnx*nny              # upper right node number
        lrn = nnx                  # lower right node number
        lln = 1                    # lower left node number
        topEdge  = np.array([range(uln, urn), range(uln+1, urn+1)]).T   # 上边边界单元
        #topEdge1 = topEdge
        botEdge  = np.array([range(lln, lrn), range(lln+1, lrn+1)]).T  # 底边边界单元
        rightEdge = np.array([range(lrn, urn-nnx+1, nnx), range(lrn+nnx, urn+1, nnx)]) 
        # GET NODES ON ESSENTIAL BOUNDARY

        botNodes = np.unique(botEdge)     # 底边边界上的节点编号
        topNodes = np.unique(topEdge)     # 上边边界上的节点编号
        rightNodes = np.unique(rightEdge) # 右边界上的节点编号
        leftNodes = rightNodes-(nnx-1)    # 左边界上的节点编号
        dispNodes = botNodes              # x和y方向均被约束
        rightNodes1 = rightNodes[1:]
        leftNodes1 = leftNodes[1:]
        dispNodes1 = np.union1d(leftNodes1, rightNodes1)  # 只有水平方向上位移被约束
    return  topEdge, dispNodes, dispNodes1, leftNodes1, topNodes

# function [topEdge,topEdge1,dispNodes,dispNodes1,leftNodes1]=supportcond(elemType,numx,numy)
# % Forms a vector of nodes which are restrained in x direction or in both
# % x and y direction which are the supports of the domain. 
# %It also forms a matrix of nodes for load application at the top of the domain.

# switch elemType
#     case 'Q9'
#         nnx=numx*2+1;
#         nny=numy*2+1;
#         urn =nnx*nny;              % upper right node number
#         uln =urn-(nnx-1);          % upper left node number
#         lrn = nnx;                 % lower right node number
#         lln = 1;                   % lower left node number

#         topEdge  = [ uln:1:(urn-1); (uln+1):1:urn ]';
#         topEdge1=[uln:2:(urn-2); (uln+1):2:(urn-1);(uln+2):2:urn ]';
#         botEdge  = [ lln:1:(lrn-1); (lln+1):1:lrn ]';
#         leftEdge=[(lln:nnx:(uln-nnx));(lln+nnx:nnx:uln)]';
#         rightEdge=[(nnx:nnx:(urn-nnx));(nnx+nnx:nnx:urn)]';
#         % GET NODES ON ESSENTIAL BOUNDARY

#         botNodes   = unique(botEdge);
#         topNodes   = unique(topEdge);
#         leftNodes = unique(leftEdge);
#         rightNodes = unique(rightEdge);
#         dispNodes = botNodes;
#         rightNodes1=rightNodes(2:end);
#         leftNodes1=leftNodes(2:end);
#         dispNodes1=union(leftNodes1,rightNodes1);
#     case 'Q8'
#         nnx=numx*2+1;
#         nny=numy*2+1;        
#         urn =(nnx*nny)-(numx*numy);% upper right node number
#         uln =urn-(nnx-1);          % upper left node number
#         lrn = nnx;                 % lower right node number
#         lln = 1;                   % lower left node number

#         topEdge  = [ uln:1:(urn-1); (uln+1):1:urn ]';
#         topEdge1=[uln:2:(urn-2); (uln+1):2:(urn-1);(uln+2):2:urn ]';
#         botEdge  = [ lln:1:(lrn-1); (lln+1):1:lrn ]';
#         % GET NODES ON ESSENTIAL BOUNDARY

#         botNodes   = unique(botEdge);
#         topNodes   = unique(topEdge);
#         leftNodes = union((lln:nnx+numx+1:uln),(nnx+1:nnx+numx+1:uln-numx+1))';
#         rightNodes = union((nnx:nnx+numx+1:urn),(nnx+numx+1:nnx+numx+1:urn-nnx))';
#         dispNodes = botNodes;
#         rightNodes1=rightNodes(2:end);
#         leftNodes1=leftNodes(2:end);
#         dispNodes1=union(leftNodes1,rightNodes1);
#         leftEdge=[leftNodes(1:nny-1,:),leftNodes1];
#         rightEdge=[rightNodes(1:nny-1,:),rightNodes1];
#     case {'Q4','T3'}
#         nnx=numx+1;
#         nny=numy+1;
#         uln = nnx*(nny-1)+1;       % upper left node number
#         urn = nnx*nny;             % upper right node number
#         lrn = nnx;                 % lower right node number
#         lln = 1;                   % lower left node number
#         topEdge  = [ uln:1:(urn-1); (uln+1):1:urn ]';
#         topEdge1=topEdge;
#         botEdge  = [ lln:1:(lrn-1); (lln+1):1:lrn ]';
#         rightEdge = (lrn:nnx:(urn))';
#         % GET NODES ON ESSENTIAL BOUNDARY

#         botNodes   = unique(botEdge);
#         topNodes   = unique(topEdge);
#         rightNodes = unique(rightEdge);
#         leftNodes = rightNodes-(nnx-1);
#         dispNodes = botNodes;
#         rightNodes1=rightNodes(2:end);
#         leftNodes1=leftNodes(2:end);
#         dispNodes1=union(leftNodes1,rightNodes1);
#     case 'T6'
#         nnx=numx*2+1;
#         nny=numy*2+1;        
#         urn =nnx*nny;              % upper right node number
#         uln =urn-(nnx-1);          % upper left node number
#         lrn = nnx;                 % lower right node number
#         lln = 1;                   % lower left node number

#         topEdge  = [ uln:1:(urn-1); (uln+1):1:urn ]';
#         topEdge1=[uln:2:(urn-2); (uln+1):2:(urn-1);(uln+2):2:urn ]';
#         botEdge  = [ lln:1:(lrn-1); (lln+1):1:lrn ]';
#         % GET NODES ON ESSENTIAL BOUNDARY

#         botNodes   = unique(botEdge);
#         topNodes   = unique(topEdge);
#         leftNodes =(lln:nnx:uln)';
#         rightNodes =(lrn:nnx:urn)';
#         dispNodes = botNodes;
#         rightNodes1=rightNodes(2:end);
#         leftNodes1=leftNodes(2:end);
#         dispNodes1=union(leftNodes1,rightNodes1);
#         leftEdge=[leftNodes(1:nny-1,:),leftNodes1];
#         rightEdge=[rightNodes(1:nny-1,:),rightNodes1];
# end

# end   % end of function
