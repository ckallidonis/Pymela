'''
Created on Dec. 11, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Class definition that performs and holds Linear fit data on the Summed ratio
'''

import pymela.io.json_io as JSONio
import pymela.io.file_formats as ioForm
import pymela.tools.tag_creators as tags
import pymela.tools.jackknife as jackknife
import pymela.tools.gamma as gmat
#import pymela.fit.linear_fit as linearFit

import numpy as np
import h5py


# The class holding the summation method fits
# The fist performed are of the form y = M*x + b, where M,b are fit parameters, and M is the desired matrix element
#
class SummationFit():
    def __init__(self, ratio, ratioType, fitInfo, analysisInfo):
        self.ratioBins = ratio.bins[ratioType]
        self.ratioMean = ratio.mean[ratioType]

        self.fitInfo = fitInfo # The list of types of fits
        self.analysisInfo = analysisInfo

        # Real-Imaginary part
        self.RI = ['Re','Im']

        # Generic definitions for each type of fit
        self.fitParams = {'Linear': ['M','b']}
    
        # What type of fits we will perform:
        self.fitPerform = []
        for fitSeq in self.fitInfo:
            if fitSeq['Type'] not in self.fitParams.keys():
                raise ValueError('Fit type %s not implemented (yet)! Supported types are: '%(fitSeq['Type']), self.fitParams.keys())
            self.fitPerform.append(fitSeq['Type'])


        # Prepare fit data scructure
        self.bins = {} # The fit parameter bins
        self.mean = {} # and mean
        self.chiBins = {} # Chi-square of the fit
        self.chiMean = {} # Chi-square of the fit

        # Each fit type has different needs and parameters, so we have to treat each one separately
        # concerning the fit data
        for fitSeq in self.fitInfo:
            fType = fitSeq['Type']
            fLabel = fitSeq['Label']
            tsepLowList = fitSeq['tsepLow']
            fPrmList = self.fitParams[fType]

            self.bins[fLabel] = {} # The fit parameter bins
            self.mean[fLabel] = {} # and mean
            self.chiBins[fLabel] = {} # Chi-square
            self.chiMean[fLabel] = {} # of the fit
            if fType == 'Linear':
                for tsepL in tsepLowList:
                    sLTag = 'tL%d'%(tsepL)
                    self.chiBins[fLabel][sLTag] = {}
                    self.chiMean[fLabel][sLTag] = {}
                    for fP in fPrmList:
                        fpTag = fP + '_%s'%(sLTag)
                        self.bins[fLabel][fpTag] = {}
                        self.mean[fLabel][fpTag] = {}
                        for ri in self.RI:
                            self.bins[fLabel][fpTag][ri] = {}
                            self.mean[fLabel][fpTag][ri] = {}
                    for ri in self.RI:
                        self.chiBins[fLabel][sLTag][ri] = {}
                        self.chiMean[fLabel][sLTag][ri] = {}
        #--------------------------

        self.momAvg = ratio.momAvg
        self.dispAvg = ratio.dispAvg
        self.Nbins = ratio.Nbins

        self.dSetAttr3pt = ratio.dSetAttr3pt
        self.dSetAttr2pt = ratio.dSetAttr2pt
    
    print('Summation Fits initialized')
    # End __init__() -------------


