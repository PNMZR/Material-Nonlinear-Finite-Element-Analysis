# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 21:03:13 2020

@author: lenovo
"""
import numpy as np

def formdg(tsi, q, theta):
    """
    calculates the partial derivatives of the plastic potential with respect
    to p, J2 and J3.
    """
    dg1 = np.sin(tsi/180*np.pi)
    if np.sin(theta/180*np.pi) > 0.49:   # close to theta=30 corner,smoothen the curve with
        sw = 1            # triaxial compression case s1 > s2 = s3
        dg2 = (0.25/q)*(3-sw*np.sin(tsi/180*np.pi))
        dg3 = 0
    
    elif -1*np.sin(theta/180*np.pi) > 0.49:   # close to theta=-30 corner,smoothen the curve
        sw = -1                 # with triaxial extension case s1 = s2 > s3
        dg2 = (0.25/q)*(3-sw*np.sin(tsi/180*np.pi))
        dg3 = 0
    else:                        # all other cases
        dg2 = (np.sqrt(3)*np.cos(theta/180*np.pi)/(2*q))*(1+(np.tan(theta/180*np.pi)*np.tan(3*theta/180*np.pi)) +
        ((np.sin(tsi/180*np.pi)/np.sqrt(3))*(np.tan(3*theta/180*np.pi)-np.tan(theta/180*np.pi))))
        dg3 = 1.5*((np.sqrt(3)*np.sin(theta/180*np.pi))+(np.sin(tsi/180*np.pi)*np.cos(theta/180*np.pi))/(q**2*np.cos(3*theta/180*np.pi)))
    return dg1, dg2, dg3
    

# function [dg1,dg2,dg3] =formdg(tsi,q,theta)

# % calculates the partial derivatives of the plastic potential with respect
# % to p, J2 and J3.

# dg1=sind(tsi);
# if sind(theta)>0.49   % close to theta=30 corner,smoothen the curve with
#     sw=1;            % triaxial compression case s1 > s2 = s3
#     dg2=(0.25/q)*(3-sw*sind(tsi));
#     dg3=0;
    
# elseif -1*sind(theta)>049  % close to theta=-30 corner,smoothen the curve
#     sw=-1;               % with triaxial extension case s1 = s2 > s3
#     dg2=(0.25/q)*(3-sw*sind(tsi));
#     dg3=0;
# else                        % all other cases
#     dg2=(sqrt(3)*cosd(theta)/(2*q))*(1+(tand(theta)*tand(3*theta))+...
#         ((sind(tsi)/sqrt(3))*(tand(3*theta)-tand(theta))));
#     dg3=1.5*((sqrt(3)*sind(theta))+(sind(tsi)*cosd(theta))/(q^2*cosd(3*theta)));
    
# end
# end % end of function