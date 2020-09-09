# -*- coding: utf-8 -*-

import numpy as np
import feglqd1

def feglqd3(nglx, ngly, nglz):
    """
    Purpose:
        Determine the integration points and weighting coefficients of
        Gauss-Legendre quadrature for three-dimensional integration.

    Synopsis: point3, weight3 = feglqd3(nglx, ngly, nglz)

    Variable Description:
        nglx - number of integration points in the x-axis
        ngly - number of integration points in the y-axis
        nglz - number of integration points in the z-axis
        point3 - vector containing integration points
        weight3 - vector containing weighting coefficients
        """

    # -----------------------------------------------------
    # determine the largest one between nglx, ngly and nglz
    # -----------------------------------------------------
    if nglx > ngly:
        if nglx > nglz:
            ngl = nglx
        else:
            ngl = nglz
    else:
        if ngly > nglz:
            ngl = ngly
        else:
            ngl = nglz

    # --------------
    # initialization
    # --------------
    point3 = np.zeros((ngl, 3))
    weight3 = np.zeros((ngl, 3))

    # -------------------------------------------------
    # find corresponding integration points and weights
    # -------------------------------------------------
    pointx, weightx = feglqd1.feglqd1(nglx)
    pointy, weighty = feglqd1.feglqd1(ngly)
    pointz, weightz = feglqd1.feglqd1(nglz)

    # ------------------------------
    # quadrature for three-dimension
    # ------------------------------
    for intx in range(nglx):
        point3[intx, 0] = pointx[intx]
        weight3[intx, 0] = weightx[intx]

    for inty in range(ngly):
        point3[inty, 1] = pointy[inty]
        weight3[inty, 1] = weighty[inty]

    for intz in range(nglz):
        point3[intz, 2] = pointz[intz]
        weight3[intz, 2] = weightz[intz]
    return point3, weight3