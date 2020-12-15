'''
Created on Dec. 15, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Class definition that computes reduced Ioffe-time distributions
'''

import pymela.io.json_io as JSONio
import pymela.io.file_formats as ioForm
import pymela.tools.tag_creators as tags
import pymela.tools.jackknife as jackknife
import pymela.tools.gamma as gmat

import numpy as np
import h5py
import scipy.optimize as scipyOpt

# The class holding the Ioffe-time Distributions
#
class ITD():
    def __init__(self, plat = None, summ = None, ITDinfo = None, fitInfo = None):

        if plat == None and summ == None:
            raise ValueError('All of the supported fit types are "None". Cannot define ITDs!')

        self.fitInfo = fitInfo
        self.info = ITDinfo

        # Real-Imaginary part
        self.RI = ['Re','Im']

        # Types of fits that we are considering for the ITDs
        self.fitLabels = self.info['Optimal Fits'].keys()

        # Make sure that the input labels are included in the fits performed earlier
        self.plat = plat
        if self.plat != None:
            for fitSeq in self.plat.fitInfo:
                if fitSeq['Label'] not in self.fitLabels:
                    raise ValueError('Fit Label %s not in Input Fit Labels'%(fitSeq['Label']))

        self.summ = summ
        if self.summ != None:
            for fitSeq in self.summ.fitInfo:
                if fitSeq['Label'] not in self.fitLabels:
                    raise ValueError('Fit Label %s not in Input Fit Labels'%(fitSeq['Label']))
        #--------------------------------------

        if self.plat != None: 
            # If self.summ also != None then that's still OK, because these attributes are the same
            # by definition in plat and summ
            self.momAvg  = self.plat.momAvg
            self.dispAvg = self.plat.dispAvg
            self.Nbins   = self.plat.Nbins
            self.dSetAttr3pt = self.plat.dSetAttr3pt
        else:
            # Get these attributes from the summ fits instead, it MUST be defined otherwise ValueError is raised
            self.momAvg  = self.summ.momAvg
            self.dispAvg = self.summ.dispAvg
            self.Nbins   = self.summ.Nbins
            self.dSetAttr3pt = self.summ.dSetAttr3pt


        # The ITD bins and mean
        self.bins = {} # The fit parameter bins
        self.mean = {} # and mean
        for fit in self.fitLabels:
            self.bins[fit] = {}
            self.mean[fit] = {}
            for ri in self.RI:
                self.bins[fit][ri] = {}
                self.mean[fit][ri] = {}
        #--------------------------
        print('ITD initialized')
    # End __init__() -------------

