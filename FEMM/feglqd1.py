# -*- coding: utf-8 -*-

import numpy as np

def feglqd1(ngl):
    """
    Purpose:
        Determine the integration points and weighting coefficients of
        Gauss-Legendre quadrature for one-dimensional integration.

    Synopsis: point1, weight1 = feglqdl(ngl)

    Variable Description:
        ngl - number of integration points
        point1 - vector containing integration points
        weight1 - vector containing weighting coefficients
    """

    # --------------
    # initialization
    # --------------
    point1 = np.zeros(ngl)
    weight1 = np.zeros(ngl)

    # -------------------------------------------------
    # find corresponding integration points and weights
    # -------------------------------------------------
    if ngl == 1:                        # 1-point quadrature rule
        point1[0] = 0.0
        weight1[0] = 2.0
    elif ngl == 2:                      # 2-point quadrature rule
        point1[0] = -0.577350269189626
        point1[1] = -point1[0]
        weight1[0] = 1.0;
        weight1[1] = weight1[0]
    elif ngl == 3:                      # 3-point quadrature rule
        point1[0] = -0.774596669241483
        point1[1]= 0.0
        point1[2] = -point1[0]
        weight1[0] = 0.555555555555556
        weight1[1] = 0.888888888888889
        weight1[2] = weight1[0]
    elif ngl == 4:                      # 4-point quadrature rule
        point1[0] = -0.861136311594053
        point1[1] = -0.339981043584856
        point1[2] = -point1[1]
        point1[3] = -point1[0]
        weight1[0] = 0.347854845137454
        weight1[1] = 0.652145154862546
        weight1[2] = weight1[1]
        weight1[3] = weight1[0]
    else:                               # 5-point quadrature rule
        point1[0] = -0.906179845938664
        point1[1] = -0.538469310105683
        point1[2] = 0.0
        point1[3] = -point1[1]
        point1[4] = -point1[0]
        weight1[0] = 0.236926885056189
        weight1[1] = 0.478628670499366
        weight1[2] = 0.568888888888889
        weight1[3] = weight1[1]
        weight1[4] = weight1[0]
    return point1, weight1