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
import pymela.fit.linear_fit as linearFit

import numpy as np
import h5py
import scipy.optimize as scipyOpt

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
        self.tsepFitX = {} # The x-axis data for each fit
        for fitSeq in self.fitInfo:
            fType = fitSeq['Type']
            fLabel = fitSeq['Label']
            tsepLowList = fitSeq['tsepLow']
            fPrmList = self.fitParams[fType]

            self.tsepFitX[fLabel] = {}

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

    def performFits(self):

        def makeLinearFit(fitSeq):
            fType = fitSeq['Type']
            fLabel = fitSeq['Label']
            tsepLowList = fitSeq['tsepLow']
            fPrmList = self.fitParams[fType]

            Nfdata = {}
            xdata  = {}

            for mom in self.momAvg:
                mTag = tags.momString(mom)
                tsepList = self.dSetAttr3pt[mTag]['tsep']
                dispListAvg = self.dispAvg[mTag]
                gammaList   = self.dSetAttr3pt[mTag]['gamma']

                # Determine the x-data for each tLow                
                self.tsepFitX[fLabel][mTag] = {}
                Nfdata[mTag] = {}
                xdata[mTag] = {}
                for tL in tsepLowList:
                    sLTag = 'tL%d'%(tL)
                    self.tsepFitX[fLabel][mTag][sLTag] = tsepList[tsepList.index(tL):]
                    xdata[mTag][sLTag] = self.tsepFitX[fLabel][mTag][sLTag]
                    Nfdata[mTag][sLTag] = len(xdata[mTag][sLTag])
                    for ri in self.RI:
                        self.chiBins[fLabel][sLTag][ri][mTag] = {}
                        self.chiMean[fLabel][sLTag][ri][mTag] = {}
                        for fP in fPrmList:
                            fpTag = fP + '_%s'%(sLTag)
                            self.bins[fLabel][fpTag][ri][mTag] = {}
                            self.mean[fLabel][fpTag][ri][mTag] = {}
        
                        for z3 in dispListAvg:
                            for gamma in gammaList:
                                dkeyF = (z3,gamma)
                                self.chiBins[fLabel][sLTag][ri][mTag][dkeyF]  = np.zeros(self.Nbins,dtype=np.float64)                                   
                                for fP in fPrmList:
                                    fpTag = fP + '_%s'%(sLTag)
                                    self.bins[fLabel][fpTag][ri][mTag][dkeyF] = np.zeros(self.Nbins,dtype=np.float64)

                                # Perform the fits for each bin
                                for b in range(self.Nbins):
                                    ydata = np.zeros(Nfdata[mTag][sLTag],dtype=np.float64)
                                    yerr  = np.zeros(Nfdata[mTag][sLTag],dtype=np.float64)

                                    # Fill in fit data
                                    for itL, tL in enumerate(self.tsepFitX[fLabel][mTag][sLTag]):
                                        dkey = (tL,z3,gamma)
                                        ydata[itL] = self.ratioBins[ri][mTag][dkey][b]
                                        yerr[itL]  = self.ratioMean[ri][mTag][dkey][1]

                                    # Perform the fit
                                    # print('xdata:', xdata[mTag][sLTag])
                                    # print('ydata:', ydata)
                                    # print('yerr:' , yerr)
                                    fprmRes, covRes = scipyOpt.curve_fit(linearFit.model, xdata[mTag][sLTag], ydata, sigma=yerr)

                                    self.chiBins[fLabel][sLTag][ri][mTag][dkeyF][b] = linearFit.chiSquare(xdata[mTag][sLTag], ydata, yerr,
                                                                                                          fprmRes[0],fprmRes[1])
                                    for ifP, fP in enumerate(fPrmList):
                                        fpTag = fP + '_%s'%(sLTag)
                                        self.bins[fLabel][fpTag][ri][mTag][dkeyF][b] = fprmRes[ifP]
                                # End for bins

                                # Jackknife averages
                                self.chiMean[fLabel][sLTag][ri][mTag][dkeyF] = jackknife.mean(self.chiBins[fLabel][sLTag][ri][mTag][dkeyF],
                                                                                              Nbins = self.Nbins, Nspl=1)
                                for fP in fPrmList:
                                    fpTag = fP + '_%s'%(sLTag)
                                    self.mean[fLabel][fpTag][ri][mTag][dkeyF] = jackknife.mean(self.bins[fLabel][fpTag][ri][mTag][dkeyF],
                                                                                                Nbins = self.Nbins, Nspl=1)
                # End for tsepLow ------
                print('%s fits for momentum %s completed'%(fType,mTag))
        # End makeLinearFit ---------

        for fitSeq in self.fitInfo:
            if fitSeq['Type'] == 'Linear':
                makeLinearFit(fitSeq)

    # End performFits() -------------

