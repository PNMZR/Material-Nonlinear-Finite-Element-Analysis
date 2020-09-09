# -*- coding: utf-8 -*-

"""
EX6.5.2.py

Gauss-Legendre quadrature of a function in 2-dimension

Problem description:
    
    Integrate f(x,y) = 1+4x*y-3x^2*y^2+x^4*y^6 over -1<x<1 and -1<y<1

Variable descriptions:
    point2 = integration (or sampling) points
    weight2 = weighting coefficients
    nglx = number of integration points along x-axis
    ngly = number of integration points along y-axis
"""

import feglqd2

nglx = 3   # (2*nglx-1)=4
ngly = 4   # (2*ngly-1)=6

point2, weight2 = feglqd2.feglqd2(nglx, ngly)   # sampling points and weights

# -----------------------------------
# summation for numerical integration
# -----------------------------------
value = 0.0
for intx in range(nglx):
    x = point2[intx, 0]                # sampling point in x-axis
    wtx = weight2[intx, 0]             # weight in x-axis
    for inty in range(ngly):
        y = point2[inty, 1]            # sampling point in y-axis
        wty = weight2[inty, 1]         # weight in y-axis
        func = 1 + 4*x*y- 3*x**2*y**2 + x**4*y**6  # evaluate function
        value = value + func*wtx*wty

print(value)