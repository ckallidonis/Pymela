'''
Created on Nov.30, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

This file contains functions related to Jackknife sampling
'''

import numpy as np

# Define the number of bins
def Nbins(Ndata, binsize=1):
    mod   = Ndata%binsize
    return (Ndata - mod) // binsize # That's always an integer
#-------------------------------------

def sampling(sample, Nbins, binsize=1):

    if len(np.shape(sample)) != 1:
        raise ValueError('Jackknife sampling: Works only with 1-d samples for now!')

    Ndata = len(sample)
    mod   = Ndata%binsize

    bins = np.zeros(Nbins,dtype=sample.dtype)

    csum = np.sum(sample) # Sum w.r.t to the configurations
    for m in np.arange(1,mod+1):
        csum -= sample[Ndata-m]  # Throw away data in case there is modulo
            
    for b in range(Nbins):
        bsum = 0
        for k in range(binsize):
            bsum += sample[b*binsize+k]
            
        bins[b] = (csum - bsum) / float(Ndata - binsize - mod) # Bin averages for each bin,t

    return bins
#-------------------------------------

def mean(bins,Nbins,Nspl):

    if np.shape(bins)[0] != Nbins:
        raise ValueError('Jackknife mean: The sampled dimension must be the first one in the "bins" array')
    Jax = 0
    
    if(Nspl==1):
        ave = np.mean(bins,dtype=bins.dtype)

        sqsum = sum(map(lambda x:x*x,ave-bins))
        fac = (Nbins -1) / float(Nbins)
        err = np.sqrt(fac*sqsum)
    else:
        ave = np.mean(bins, axis=Jax,dtype=bins.dtype)

        err = np.zeros(Nspl,dtype=np.float128)
        for i in range(Nspl):
            sqsum = sum(map(lambda x:x*x,ave[i]-bins[:,i]))
            fac = (Nbins -1) / float(Nbins)
            err[i] = np.sqrt(fac*sqsum.real)        
            
    return (ave,err)
#-------------------------------------
