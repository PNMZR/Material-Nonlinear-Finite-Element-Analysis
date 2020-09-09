# -*- coding: utf-8 -*-

import numpy as np

def feode2l(acoef, bcoef, ccoef, eleng):
    """
    Purpose:
    Element matrix for (a u" + b u'+ c u) using linear element

    Synopsis:
    k = feode2l(acoef, bcoef, ccoef, eleng)

    Variable Description:
    k - element matrix (size of 2x2)
    acoef - coefficient of the second order derivative term
    bcoef - coefficient of the first order derivative term
    ccoef - coefficient of the zero-th order derivative term
    """

    al = -(acoef/eleng)
    a2 = bcoef/2
    a3 = ccoef*eleng/6
    k = np.array([[al-a2+2*a3, -al+a2+a3], [-al-a2+a3, al+a2+2*a3]])
    return k