'''
Created on Dec. 10, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Class definition that performs and holds Plateau (Constant) fit data
'''

import pymela.io.json_io as JSONio
import pymela.io.file_formats as ioForm
import pymela.tools.tag_creators as tags
import pymela.tools.jackknife as jackknife
import pymela.tools.gamma as gmat
import pymela.fit.constant_fit as constFit

import numpy as np
import h5py


# The class holding the plateau fits
#
class PlateauFit():
    def __init__(self, ratio, ratioType, fitInfo, analysisInfo):
        self.ratioBins = ratio.bins[ratioType]
        self.ratioMean = ratio.mean[ratioType]

        self.fitInfo = fitInfo
        self.analysisInfo = analysisInfo

        # Real-Imaginary part
        self.RI = ['Re','Im']

        self.Mbins = {} # The matrix element (constant fit value)
        self.Mmean = {} # The matrix element (constant fit value)
        self.chiBins = {} # Chi-square of the fit
        self.chiMean = {} # Chi-square of the fit
        self.optimalFit = {} # Structure that holds the optimal plateau fits
        for fitSeq in self.fitInfo:
            fType = fitSeq['Type']
            fLabel = fitSeq['Label']
            if fType != 'Constant':
                print('PlateauFits: Supports only Constant fits for now!')
            self.Mbins[fLabel] = {}
            self.Mmean[fLabel] = {}
            self.chiBins[fLabel] = {}
            self.chiMean[fLabel] = {}
            self.optimalFit[fLabel] = {}
            for ri in self.RI:
                self.Mbins[fLabel][ri] = {}
                self.Mmean[fLabel][ri] = {}
                self.chiBins[fLabel][ri] = {}
                self.chiMean[fLabel][ri] = {}
                self.optimalFit[fLabel][ri] = {}

        self.momAvg = ratio.momAvg
        self.dispAvg = ratio.dispAvg
        self.Nbins = ratio.Nbins

        self.dSetAttr3pt = ratio.dSetAttr3pt
        self.dSetAttr2pt = ratio.dSetAttr2pt

        # Define required fit structures
        self.fitAttr = {}
        for mom in self.momAvg:
            mTag = tags.momString(mom)
            tsepList = self.dSetAttr3pt[mTag]['tsep']

            self.fitAttr[mTag] = {}
            for tsep in tsepList:
                self.fitAttr[mTag][tsep] = {}

                self.fitAttr[mTag][tsep]['Nfits'] = tsep//2 - 1 # How many fits for each tsep

                tini,tfin = 1, tsep-1 # Full range, omit the source point, go up to the end
                self.fitAttr[mTag][tsep]['Rng'] = []
                for nf in range(self.fitAttr[mTag][tsep]['Nfits']):
                    self.fitAttr[mTag][tsep]['nf=%d'%(nf)] = {}
                    tstart, tstop = tini+nf,tfin-nf  # Range of each fit
                    Npts = tstop-tstart+1  # How many points in each fit
                    self.fitAttr[mTag][tsep]['nf=%d'%(nf)]['xdata']  = np.arange(Npts+1)
                    self.fitAttr[mTag][tsep]['nf=%d'%(nf)]['tstart'] = tstart # Range
                    self.fitAttr[mTag][tsep]['nf=%d'%(nf)]['tstop']  = tstop  # of each fit
                    self.fitAttr[mTag][tsep]['nf=%d'%(nf)]['Npts']   = Npts   # How many points in each fit
                    self.fitAttr[mTag][tsep]['Rng'].append('%d-%d'%(tstart,tstop))                    
    
    print('Plateau Fits initialized')
    # End __init__() -------------

    def performFits(self):

        def makeConstantFits(fitSeq):
            fType = fitSeq['Type']
            fLabel = fitSeq['Label']
            chiCrit = fitSeq['Chi Criterion']

            for mom in self.momAvg:
                mTag = tags.momString(mom)
                tsepList = self.dSetAttr3pt[mTag]['tsep']
                dispListAvg = self.dispAvg[mTag]
                gammaList   = self.dSetAttr3pt[mTag]['gamma']
                
                for ri in self.RI:
                    self.Mbins[fLabel][ri][mTag] = {}
                    self.Mmean[fLabel][ri][mTag] = {}
                    self.chiBins[fLabel][ri][mTag] = {}
                    self.chiMean[fLabel][ri][mTag] = {}
                    self.optimalFit[fLabel][ri][mTag] = {}

                    for tsep in tsepList:
                        fAttr = self.fitAttr[mTag][tsep]
                        for z3 in dispListAvg:
                            for gamma in gammaList:
                                dkey = (tsep,z3,gamma)

                                self.Mbins[fLabel][ri][mTag][dkey] = {}
                                self.Mmean[fLabel][ri][mTag][dkey] = {}
                                self.chiBins[fLabel][ri][mTag][dkey] = {}
                                self.chiMean[fLabel][ri][mTag][dkey] = {}
                                self.optimalFit[fLabel][ri][mTag][dkey] = -1 # Negative number means no Optimal Fit Found
                                
                                optimalFitFound = False
                                for nf in range(fAttr['Nfits']):
                                    tstart = fAttr['nf=%d'%(nf)]['tstart']
                                    tstop  = fAttr['nf=%d'%(nf)]['tstop']

                                    self.Mbins[fLabel][ri][mTag][dkey][nf] = np.zeros(self.Nbins,dtype=np.float128)
                                    self.chiBins[fLabel][ri][mTag][dkey][nf] = np.zeros(self.Nbins,dtype=np.float128)

                                    for b in range(self.Nbins):
                                        data = self.ratioBins[ri][mTag][dkey][b,tstart:tstop+1]
                                        err  = self.ratioMean[ri][mTag][dkey][1][tstart:tstop+1]
                                        if fType == 'Constant':
                                            self.Mbins[fLabel][ri][mTag][dkey][nf][b]   = constFit.fit(data,err)
                                            self.chiBins[fLabel][ri][mTag][dkey][nf][b] = constFit.chiSquare(data,err,self.Mbins[fLabel][ri][mTag][dkey][nf][b])

                                    self.Mmean[fLabel][ri][mTag][dkey][nf]   = jackknife.mean(self.Mbins[fLabel][ri][mTag][dkey][nf],
                                                                                              Nbins = self.Nbins, Nspl=1)
                                    self.chiMean[fLabel][ri][mTag][dkey][nf] = jackknife.mean(self.chiBins[fLabel][ri][mTag][dkey][nf],
                                                                                              Nbins = self.Nbins, Nspl=1)

                                    # Determine optimal plateau fit
                                    if (self.chiMean[fLabel][ri][mTag][dkey][nf][0] <= chiCrit) and not optimalFitFound:
                                        self.optimalFit[fLabel][ri][mTag][dkey] = nf
                                        optimalFitFound = True
                                # End for Nfits
                print('%s fits, with label %s for momentum %s completed.'%(fType, fLabel, mTag))
        # End makeFits() ----------------

        for fitSeq in self.fitInfo:
            if fitSeq['Type'] == 'Constant':
                makeConstantFits(fitSeq)
    # End performFits() -------------

    def writeHDF5(self):

        def dumpHDF5(fitSeq,h5_file):
            fType = fitSeq['Type']
            fLabel = fitSeq['Label']

            for mom in self.momAvg:
                mTag = tags.momString(mom)
                mh5Tag = tags.momH5(mom)
                tsepList = self.dSetAttr3pt[mTag]['tsep']
                dispListAvg = self.dispAvg[mTag]
                gammaList   = self.dSetAttr3pt[mTag]['gamma']
                
                for z3 in dispListAvg:
                    dispTag = tags.disp(z3)
                    for gamma in gammaList:
                        insTag = tags.insertion(gamma)
                        for tsep in tsepList:
                            dkey = (tsep,z3,gamma)
                            tsepTag = tags.tsep(tsep)
                            fAttr = self.fitAttr[mTag][tsep]

                            for ri in self.RI:
                                # Write optimalFitValues
                                group = '%s/%s/%s/%s/%s'%(ri,mh5Tag,dispTag,insTag,tsepTag)
                                dset_name = 'OptimalFitRanges/' + group
                                h5_file.create_dataset(dset_name, data = np.array([self.optimalFit[fLabel][ri][mTag][dkey]]))

                                for nf in range(fAttr['Nfits']):
                                    tstart = fAttr['nf=%d'%(nf)]['tstart']
                                    tstop  = fAttr['nf=%d'%(nf)]['tstop']
                                    h5LabelNf = 'nf%d_%d-%d'%(nf,tstart,tstop)

                                    group = '%s/%s/%s/%s/%s/%s'%(ri,mh5Tag,dispTag,insTag,tsepTag,h5LabelNf)
                                    dset_name_Mbins = 'MatElem/bins/' + group 
                                    dset_name_Mmean = 'MatElem/mean/' + group 
                                    dset_name_chiBins = 'chiSquare/bins/' + group 
                                    dset_name_chiMean = 'chiSquare/mean/' + group 

                                    h5_file.create_dataset(dset_name_Mbins, data = self.Mbins[fLabel][ri][mTag][dkey][nf])
                                    h5_file.create_dataset(dset_name_Mmean, data = self.Mmean[fLabel][ri][mTag][dkey][nf],dtype='f')
                                    h5_file.create_dataset(dset_name_chiBins, data = self.chiBins[fLabel][ri][mTag][dkey][nf])
                                    h5_file.create_dataset(dset_name_chiMean, data = self.chiMean[fLabel][ri][mTag][dkey][nf],dtype='f')
            # End for momentum
            print('Plateau fitting data for type = %s, label = %s written in HDF5.'%(fType,fLabel))
        # End dumpHDF5 ----------------

        for fitSeq in self.fitInfo:
            if fitSeq['Write HDF5 Output']:
                h5_file = h5py.File(fitSeq['HDF5 Output File'],'w')
                dumpHDF5(fitSeq,h5_file)
                h5_file.close()
    # End writeHDF5() -------------
