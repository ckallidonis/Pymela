'''
Created on Dec.4, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Module that contains functions related to Constant fits
'''

import numpy as np

def fit(data,err):
    Sy = sum( map(lambda x:x,data/(err**2)) )    # Sy = SUM_i data[i]/err[i]**2
    S  = sum( map(lambda x:1.0/(x*x),err) )      # S  = SUM_i ( 1/err[i] )^2
    return Sy/S


def chiSquare(data,err,fVal):
    Ndof = np.shape(data)[0] - 2 # Degrees of freedom = Ndata - Nfit_param - 1
    return sum( map(lambda t:t*t, (data-fVal)/err) ) / Ndof
