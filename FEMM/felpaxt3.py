# -*- coding: utf-8 -*-

import numpy as np

def felpaxt3(r, z):
    """
    Purpose:
    Element matrix for axisymmetric Laplace equation using three-node linear
    triangular element

    Synopsis: k = felpaxt3(r, z)

    Variable Description:
        k - element stiffness matrix (size of 3x3)
        r, z - r and z coordinate values in the present element
    """
    
    k = np.zeros((3, 3))
    r1 = r[0]
    r2 = r[1]
    r3 = r[2]
    z1 = z[0]
    z2 = z[1]
    z3 = z[2]
    # area of the triangle
    A = 0.5*(r2*z3 + r1*z2 + r3*z1 - r2*z1 - r1*z3- r3*z2) 
    rc = (r1+r2+r3) / 3  # r coordinate value of the element centroid
    k[0, 0] = ((r3-r2)*(r3-r2) + (z2-z3)*(z2-z3)) / (4*A)
    k[0, 1] = ((r3-r2)*(r1-r3) + (z2-z3)*(z3-z1)) / (4*A)
    k[0, 2] = ((r3-r2)*(r2-r1) + (z2-z3)*(z1-z2)) / (4*A)
    k[1, 0] = k[0, 1]
    k[1, 1] = ((r1-r3)*(r1-r3) + (z3-z1)*(z3-z1)) / (4*A)
    k[1, 2] = ((r1-r3)*(r2-r1) + (z3-z1)*(z1-z2)) / (4*A)
    k[2, 0] = k[0, 2]
    k[2, 1] = k[1, 2]
    k[2, 2] = ((r2-r1)*(r2-r1) + (z1-z2)*(z1-z2)) / (4*A)
    k = 2*np.pi*rc*k
    return k