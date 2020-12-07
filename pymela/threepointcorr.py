'''
Created on Nov.20, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Class definition that holds data and manipulates Three-point correlation functions
'''

import pymela.io.json_io as JSONio
import pymela.io.file_formats as ioForm
import pymela.tools.tag_creators as tags
import pymela.tools.jackknife as jackknife
import pymela.tools.gamma as gmat


import numpy as np
import h5py


# The class holding the three-point correlation functions
#
# Three-point functions have dependence on the following attributes:
#   the Momentum
#   the source time, t0
#   the source and sink operators
#   the rows of the source-sink operators
#   the insertion operator
#   the source-sink time separation
#   the displacement z3
#
# For each value of the momentum, the software will average over all t0s, source/sink operators and rows and displacements.
#
class ThreePointCorrelator():
    def __init__(self, dataInfo, analysisInfo):
        self.dataInfo = dataInfo
        self.analysisInfo = analysisInfo

        #   The three-point function also has a significant real and imaginary part
        self.RI = ['Re','Im']

        # Data containers        
        # The data that is read/loaded
        # "plain" means not averaged over t0,src-snk operators, rows, or momentum, i.e. there's depedence on these attributes
        self.plainData = {}   
        self.plainBins = {}
        self.plainMean = {}

        # The averaged data over t0, src-snk operators and rows
        self.avgData = {}
        self.avgBins = {}
        self.avgMean = {}

        # The momentum- and z3-averaged data, that will be used throughout the analysis
        self.data = {}
        self.bins = {}
        self.mean = {}

        for ri in self.RI:
            self.plainData[ri] = {}
            self.plainBins[ri] = {}
            self.plainMean[ri] = {}
            self.avgData[ri] = {}
            self.avgBins[ri] = {}
            self.avgMean[ri] = {}
            self.data[ri] = {}
            self.bins[ri] = {}
            self.mean[ri] = {}

        # The number of Jackknife bins, and the binsize (same for plain and averaged data)
        self.Nbins = 0     
        self.binsize = self.analysisInfo['Binsize']

        self.dataLoaded = False

        # Fill in Attributes
        self.Nvec = self.analysisInfo['Nvec']
        self.phaseTag = self.analysisInfo['Phasing Tag']

        self.moms  = []
        self.dSetAttr = {}
        self.dSetList = self.dataInfo['Datasets']
        for dSet in self.dSetList:
            momVec = dSet['mom']
            mTag = tags.momString(momVec) # Dataset Attributes are listed for each momentum
            self.dSetAttr[mTag] = {}

            self.moms.append(momVec)

            if dSet['Compute X-rows']:
                raise ValueError('Does not support doing cross-rows in two-point function for now!')

            for attr in ['t0','Ncfg','tsep','disp','Nrows','Compute X-rows']:
                self.dSetAttr[mTag][attr] = dSet[attr]
            self.dSetAttr[mTag]['gamma'] = dSet['Insertion Operators']


            # Read source-sink operators
            intOpFile = dSet['Interpolating Operators File']            
            self.dSetAttr[mTag]['intOpList'] = []
            with open(intOpFile) as fp:
                ops = fp.readlines()
                for op in ops:
                    self.dSetAttr[mTag]['intOpList'].append((op.split()[0],op.split()[1]))
            self.dSetAttr[mTag]['Nop'] = len(ops)

        # Get the momenta that will be averaged over
        self.momAvg = self.dataInfo['Momentum Average']

    # End __init__() -------------

    def printInfo(self):
        JSONio.dumpDictObject(self.dataInfo, '\nThreePointCorrelator - Got the following 3pt data Info:')

        print('\nParsed the following momenta:')
        print(self.moms)

        JSONio.dumpDictObject(self.dSetAttr, '\nThreePointCorrelator - Parsed the following Attributes:')        
    #-------------------------------

    def getData(self):
        dataSource = self.dataInfo['Data Source']
        
        if dataSource == 'ASCII':
            print('\nWill read data from ASCII files')

            for mom in self.moms:
                mTag = tags.momString(mom)
                mFTag = tags.momFile(mom)
                t0List = self.dSetAttr[mTag]['t0']
                tsepList = self.dSetAttr[mTag]['tsep']
                dispList = self.dSetAttr[mTag]['disp']
                Nrows = self.dSetAttr[mTag]['Nrows']
                Ncfg = self.dSetAttr[mTag]['Ncfg']

                # Determine the Jackknife sampling number of Bins
                self.Nbins = jackknife.Nbins(Ncfg,self.binsize)

                for ri in self.RI:
                    self.plainData[ri][mTag] = {}

                for tsep in tsepList:
                    Nt = tsep
                    tsepTag = tags.tsep(tsep)
                    for t0 in t0List:
                        t0Tag = tags.t0(t0)
                        fileDir = ioForm.getThreePointDirASCII(self.dataInfo['Data Main Directory'],t0Tag,tsepTag,mFTag)
                        print('Reading three-point data for momentum %s, tsep = %d, t0 = %d'%(mTag,tsep,t0))
                        for z3 in dispList:
                            dispTag = tags.disp(z3)
                            for iop,opPair in enumerate(self.dSetAttr[mTag]['intOpList']):
                                srcOp,snkOp = opPair
                                for row in range(1,Nrows+1):

                                    for gamma in self.dSetAttr[mTag]['gamma']:
                                        dkey = (tsep,t0,z3,iop,row,gamma)

                                        # Determine gamma matrix name and row
                                        insOp,insRow = gmat.insertionMap(gamma)

                                        fileName = ioForm.getThreePointFileNameASCII(self.phaseTag,t0Tag,tsepTag,
                                                                                     srcOp,snkOp,row,insOp,insRow,
                                                                                     mFTag,dispTag,self.Nvec)
                                        fileRead = '%s/%s'%(fileDir,fileName)
                                        
                                        rawData = np.zeros((Ncfg,Nt), dtype=np.complex128)
                                        with open(fileRead) as fp:
                                            line = fp.readlines()
                                            c = 0
                                            for n in line:
                                                it   = c%Nt
                                                icfg = c//Nt
                                                rawData[icfg,it] = complex(np.float64(n.split()[1]),
                                                                           np.float64(n.split()[2]))
                                                c += 1
                                        # Done reading file

                                        self.plainData['Re'][mTag][dkey] = rawData.real
                                        self.plainData['Im'][mTag][dkey] = rawData.imag

                print('Reading two-point data for momentum %s completed.\n'%(mTag))

            self.dataLoaded = True
        else:
            raise ValueError('\nUnsupported "Data Source" = %s ' % (dataSource))
    # End getData() -------------


    def doStatistics(self):

        if not self.dataLoaded:
            raise ValueError('Data must be loaded first, before doing Statistical Sampling')

        for mom in self.moms:
            mTag = tags.momString(mom)
            t0List = self.dSetAttr[mTag]['t0']
            tsepList = self.dSetAttr[mTag]['tsep']
            dispList = self.dSetAttr[mTag]['disp']
            Nrows = self.dSetAttr[mTag]['Nrows']
            Ncfg = self.dSetAttr[mTag]['Ncfg']
            Nt0 = len(t0List)
            Nop = self.dSetAttr[mTag]['Nop']

            Navg = Nrows * Nt0 * Nop

            # The plain data Bins and Mean
            for ri in self.RI:
                self.plainBins[ri][mTag] = {}
                self.plainMean[ri][mTag] = {}
                self.avgData[ri][mTag] = {}
                self.avgBins[ri][mTag] = {}
                self.avgMean[ri][mTag] = {}

            for tsep in tsepList:
                Nt = tsep
                for z3 in dispList:
                    for gamma in self.dSetAttr[mTag]['gamma']:
                        dkeyAvg = (tsep,z3,gamma)
                        for ri in self.RI:
                            self.avgData[ri][mTag][dkeyAvg] = np.zeros((Ncfg,Nt),dtype=np.float64)

                        for t0 in t0List:
                            for iop,opPair in enumerate(self.dSetAttr[mTag]['intOpList']):
                                for row in range(1,Nrows+1):
                                    dkey = (tsep,t0,z3,iop,row,gamma)

                                    for ri in self.RI:
                                        # Jackknife sampling on the Plain data
                                        self.plainBins[ri][mTag][dkey] = np.zeros((self.Nbins,Nt), dtype=np.float64)
                                        for t in range(Nt):
                                            self.plainBins[ri][mTag][dkey][:,t] = jackknife.sampling(self.plainData[ri][mTag][dkey][:,t], self.Nbins, self.binsize)
                                        self.plainMean[ri][mTag][dkey] = jackknife.mean(self.plainBins[ri][mTag][dkey], self.Nbins, Nspl=Nt)

                                        # Average over Source-Sink operators, t0's and rows
                                        self.avgData[ri][mTag][dkeyAvg] += self.plainData[ri][mTag][dkey]

                        # Average over Source-Sink operators, t0's and rows
                        for ri in self.RI:
                            self.avgData[ri][mTag][dkeyAvg] = self.avgData[ri][mTag][dkeyAvg] / Navg

                            # Jackknife sampling over the averaged data, for each momentum, tsep, z3 and gamma
                            self.avgBins[ri][mTag][dkeyAvg] = np.zeros((self.Nbins,Nt), dtype=np.float64)
                            for t in range(Nt):
                                self.avgBins[ri][mTag][dkeyAvg][:,t] = jackknife.sampling(self.avgData[ri][mTag][dkeyAvg][:,t],
                                                                                          self.Nbins, self.binsize)
                            self.avgMean[ri][mTag][dkeyAvg] = jackknife.mean(self.avgBins[ri][mTag][dkeyAvg],
                                                                             self.Nbins, Nspl=Nt)

            print('\nJackknife analysis for momentum %s completed'%(mTag))
