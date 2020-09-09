# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 20:20:40 2020

@author: lenovo
"""
import numpy as np

def invariants2(kk,iel,stress):
    """
    calculates the stress invariants at each gauss point.This can be used 
    inside the iteration loop.
    """
    s_xx = stress[0,kk-1,iel-1]
    s_yy = stress[1,kk-1,iel-1]
    t_xy = stress[2,kk-1,iel-1]
    s_zz = stress[3,kk-1,iel-1]
    p = (s_xx+s_yy+s_zz)/3
    t = np.sqrt((s_xx-s_yy)**2+(s_yy-s_zz)**2+(s_zz-s_xx)**2+6*t_xy**2)/np.sqrt(3)
    q = np.sqrt(1.5)*t
    sx = s_xx - p 
    sy = s_yy - p
    sz = s_zz - p
    if q < 1e-6:
        theta = 0
    else:
        j3 = sx*sy*sz-(sx*t_xy**2)
        sine =-3*j3*np.sqrt(6)/t**3
        if sine > 1:
            sine = 1
        if sine < -1:
            sine = -1
        theta = 1/3*(np.arcsin(sine)/np.pi*180)
    return p, q, theta

# function [p,q,theta] =invariants2(kk,iel,stress)

# % calculates the stress invariants at each gauss point.This can be used  
# % inside the iteration loop.

#   s_xx=stress(1,kk,iel);s_yy=stress(2,kk,iel);
#   t_xy=stress(3,kk,iel);s_zz=stress(4,kk,iel);
#   p=(s_xx+s_yy+s_zz)/3;
#   t=sqrt((s_xx-s_yy)^2+(s_yy-s_zz)^2+(s_zz-s_xx)^2+6*t_xy^2)/sqrt(3);
#   q=sqrt(1.5)*t;
#   sx=s_xx-p; 
#   sy=s_yy-p;
#   sz=s_zz-p;
#       if q<1e-6
#          theta=0;
#      else
#          j3=sx*sy*sz-(sx*t_xy^2);
#          sine=-3*j3*sqrt(6)/t^3;
#                 if sine>1;
#                    sine=1;
#                 end
#                 if sine<-1;
#                    sine=-1;
#                 end
#           theta=1/3*(asind(sine));
#       end            

# end   % end of function