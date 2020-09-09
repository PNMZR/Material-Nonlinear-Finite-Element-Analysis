# -*- coding: utf-8 -*-

import numpy as np
import feglqd1

def feglqd2(nglx, ngly):
    """
    Purpose:
        Determine the integration points and weighting coefficients of
        Gauss-Legendre quadrature for two-dimensional integration.
    
    Synopsis: point2, weight2 = feglqd2(nglx, ngly)

    Variable Description:
        nglx - number of integration points in the x-axis
        ngly - number of integration points in the y-axis
        point2 - vector containing integration points
        weight2 - vector containing weighting coefficients
    """

    # -----------------------------------------------
    # determine the largest one between nglx and ngly
    # -----------------------------------------------
    if nglx > ngly:
        ngl = nglx
    else:
        ngl = ngly

    # --------------
    # initialization
    # --------------
    point2 = np.zeros((ngl, 2))
    weight2 = np.zeros((ngl, 2))

    # -------------------------------------------------
    # find corresponding integration points and weights
    # -------------------------------------------------
    pointx, weightx = feglqd1.feglqd1(nglx)
    pointy, weighty = feglqd1.feglqd1(ngly)

# quadrature for two-dimension
    for intx in range(nglx):
        point2[intx, 0] = pointx[intx]
        weight2[intx, 0] = weightx[intx]
    for inty in range(ngly):
        point2[inty, 1] = pointy[inty]
        weight2[inty, 1] = weighty[inty]
    return  point2, weight2