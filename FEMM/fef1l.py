# -*- coding: utf-8 -*-

import numpy as np

def fef1l(eleng):
    """
    Purpose:
    Element vector for f(x)=l using linear element
    
    Synopsis: f = fefll(xl,xr)

    Variable Description:
        f - element vector (size of 2x1)
        eleng - element length
    """
    f = np.array([eleng/2, eleng/2])
    return f