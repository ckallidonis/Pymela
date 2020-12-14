'''
Created on Dec.11, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Module that contains functions related to Linear fits
'''

import numpy as np

# The model of a linear function, y = b + M*x
def model(x,M,b):
    return b + M*x

# Chi-square for Linear fit
def chiSquare(x, y, err, M, b):
    if np.shape(x) != np.shape(y) != np.shape(err):
        raise ValueError('linear - chiSquare: Got inconsistent data shapes')
    Ndata = np.shape(x)[0]
    Ndof = Ndata - 3 # Degrees of freedom = Ndata - Nfit_param - 1
    chi = 0
    for i in range(Ndata):
        chi += ( (y[i] - model(x[i],M,b))/err[i] )**2
    return chi/Ndof
#    return sum( map(lambda t:t*t, (y[i]-model(x,M,b))/err) ) / Ndof
#------------------------
