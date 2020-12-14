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
        self.fitParams   = {'Linear': ['M','b']}
        self.fitParamsH5 = {'Linear': ['MatElem','Intersection']}
    
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

        # The structure that holds the fit bands
        self.fitBands = {}

        # Each fit type has different needs and parameters, so we have to treat each one separately
        # concerning the fit data
        self.tsepFitX = {} # The x-axis data for each fit
        for fitSeq in self.fitInfo:
            fType = fitSeq['Type']
            fLabel = fitSeq['Label']
            tsepLowList = fitSeq['tsepLow']
            fPrmList = self.fitParams[fType]

            self.tsepFitX[fLabel] = {}

            self.bins[fLabel] = {}
            self.mean[fLabel] = {}
            self.chiBins[fLabel] = {}
            self.chiMean[fLabel] = {}
            self.fitBands[fLabel] = {}

            if fType == 'Linear':
                for tsepL in tsepLowList:
                    sLTag = 'tL%d'%(tsepL)
                    self.chiBins[fLabel][sLTag] = {}
                    self.chiMean[fLabel][sLTag] = {}
                    self.fitBands[fLabel][sLTag] = {}
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
                        self.fitBands[fLabel][sLTag][ri] = {}
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

            for mom in self.momAvg:
                mTag = tags.momString(mom)
                tsepList = self.dSetAttr3pt[mTag]['tsep']
                dispListAvg = self.dispAvg[mTag]
                gammaList   = self.dSetAttr3pt[mTag]['gamma']

                # Determine the x-data for each tLow                
                self.tsepFitX[fLabel][mTag] = {}
                for tL in tsepLowList:
                    sLTag = 'tL%d'%(tL)
                    self.tsepFitX[fLabel][mTag][sLTag] = tsepList[tsepList.index(tL):]
                    xData = self.tsepFitX[fLabel][mTag][sLTag]
                    Nfdata = len(xData)
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
                                    ydata = np.zeros(Nfdata,dtype=np.float64)
                                    yerr  = np.zeros(Nfdata,dtype=np.float64)

                                    # Fill in fit data
                                    for itL, tL in enumerate(self.tsepFitX[fLabel][mTag][sLTag]):
                                        dkey = (tL,z3,gamma)
                                        ydata[itL] = self.ratioBins[ri][mTag][dkey][b]
                                        yerr[itL]  = self.ratioMean[ri][mTag][dkey][1]

                                    # Perform the fit
                                    fprmRes, covRes = scipyOpt.curve_fit(linearFit.model, xData, ydata, sigma=yerr)

                                    self.chiBins[fLabel][sLTag][ri][mTag][dkeyF][b] = linearFit.chiSquare(xData, ydata, yerr,
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


    def constructFitBands(self):

        def makeLinearFitBands(fitSeq):
            Npts = fitSeq['Fit Bands']['Npoints']

            fType = fitSeq['Type']
            fLabel = fitSeq['Label']
            tsepLowList = fitSeq['tsepLow']

            for mom in self.momAvg:
                mTag = tags.momString(mom)
                dispListAvg = self.dispAvg[mTag]
                gammaList   = self.dSetAttr3pt[mTag]['gamma']

                for tL in tsepLowList:
                    sLTag = 'tL%d'%(tL)           
                    MTag = 'M' + '_%s'%(sLTag)
                    bTag = 'b' + '_%s'%(sLTag)

                    # Determine beginning and ending of bands, and interval
                    xStart = self.tsepFitX[fLabel][mTag][sLTag][0]-1
                    xEnd   = self.tsepFitX[fLabel][mTag][sLTag][-1]+1
                    dx = (xEnd-xStart) / (Npts-1)

                    for ri in self.RI:
                        self.fitBands[fLabel][sLTag][ri][mTag] = {}
                        for z3 in dispListAvg:
                            for gamma in gammaList:
                                dkeyF = (z3,gamma)

                                self.fitBands[fLabel][sLTag][ri][mTag][dkeyF] = {'x': np.zeros(Npts,dtype=np.float64), # x
                                                                                 'v': np.zeros(Npts,dtype=np.float64), # value
                                                                                 'e': np.zeros(Npts,dtype=np.float64)} # error
                                for ix in range(Npts):
                                    x = xStart + ix*dx # Current point in the band
                                    Mmean = self.mean[fLabel][MTag][ri][mTag][dkeyF][0] # Matrix element (slope)
                                    bmean = self.mean[fLabel][bTag][ri][mTag][dkeyF][0] # Intersection

                                    self.fitBands[fLabel][sLTag][ri][mTag][dkeyF]['x'] = x
                                    self.fitBands[fLabel][sLTag][ri][mTag][dkeyF]['v'] = linearFit.model(x,Mmean,bmean)

                                    # Determine error band
                                    errBand = np.zeros(self.Nbins,dtype=np.float64)
                                    for ib in range(self.Nbins):
                                        Mbins = self.bins[fLabel][MTag][ri][mTag][dkeyF][ib]
                                        bbins = self.bins[fLabel][bTag][ri][mTag][dkeyF][ib]
                                        errBand[ib] = linearFit.model(x,Mbins,bbins)
                                    self.fitBands[fLabel][sLTag][ri][mTag][dkeyF]['e'] = jackknife.mean(errBand,self.Nbins,Nspl=1)[1]
                # End for tsepLow ------
                print('%s error bands for momentum %s completed'%(fType,mTag))
        # End makeLinearFitBands() -------------

        for fitSeq in self.fitInfo:
            if fitSeq['Type'] == 'Linear' and fitSeq['Fit Bands']['Evaluate']:
                makeLinearFitBands(fitSeq)

    # End constructFitBands() -------------

    def writeHDF5(self):

        def dumpLinearFitsHDF5(fitSeq,h5_file):
            fType = fitSeq['Type']
            fLabel = fitSeq['Label']
            tsepLowList = fitSeq['tsepLow']

            for mom in self.momAvg:
                mTag = tags.momString(mom)
                mh5Tag = tags.momH5(mom)
                dispListAvg = self.dispAvg[mTag]
                gammaList   = self.dSetAttr3pt[mTag]['gamma']
                
                for z3 in dispListAvg:
                    dispTag = tags.disp(z3)
                    for gamma in gammaList:
                        insTag = tags.insertion(gamma)
                        dkeyF = (z3,gamma)

                        for ri in self.RI:
                            for tL in tsepLowList:
                                sLTag = 'tL%d'%(tL)                                
                                tini = self.tsepFitX[fLabel][mTag][sLTag][0]
                                tfin = self.tsepFitX[fLabel][mTag][sLTag][-1]
                                h5LabelT = 'tsep_%d-%d'%(tini,tfin)

                                # Write Chi^2
                                group = '%s/%s/%s/%s/%s'%(ri,mh5Tag,dispTag,insTag,h5LabelT)
                                dset_name_chiBins = 'chiSquare/bins/' + group 
                                dset_name_chiMean = 'chiSquare/mean/' + group
                                h5_file.create_dataset(dset_name_chiBins, data = self.chiBins[fLabel][sLTag][ri][mTag][dkeyF])
                                h5_file.create_dataset(dset_name_chiMean, data = self.chiMean[fLabel][sLTag][ri][mTag][dkeyF],dtype='f')

                                # Write Fit parameters
                                for fP,fpH5 in zip(self.fitParams[fType],self.fitParamsH5[fType]):
                                    fpTag = fP + '_%s'%(sLTag)

                                    dset_name_bins = '%s/bins/'%(fpH5) + group 
                                    dset_name_mean = '%s/mean/'%(fpH5) + group
                                    h5_file.create_dataset(dset_name_bins, data = self.bins[fLabel][fpTag][ri][mTag][dkeyF])
                                    h5_file.create_dataset(dset_name_mean, data = self.mean[fLabel][fpTag][ri][mTag][dkeyF],dtype='f')
            # End for momentum
            print('Summation fitting data for type = %s, label = %s written in HDF5.'%(fType,fLabel))
        # End dumpLinearFitsHDF5 ----------------

        for fitSeq in self.fitInfo:
            if fitSeq['Write HDF5 Output']:
                h5_file = h5py.File(fitSeq['HDF5 Output File'],'w')
                if fitSeq['Type'] == 'Linear':
                    dumpLinearFitsHDF5(fitSeq,h5_file)
                h5_file.close()
    # End writeHDF5() -------------
