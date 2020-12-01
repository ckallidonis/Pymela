'''
Created on Nov.20, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Class definition that holds data and manipulates Two-point correlation functions
'''

import io.json_io as JSONio
import io.file_formats as ioForm
import tools.tag_creators as tags
import tools.jackknife as jackknife

import numpy as np

def TwoPointKeyGen(mom,t0=None,isrcOp=None,isnkOp=None,iRow=None):
    mTag = tags.momentum(mom)
    if t0==None and isrcOp==None and isnkOp==None and iRow==None:
        return (mTag)
    else:
        return (mTag,t0,isrcOp,isnkOp,iRow)
#-------------------------------

    
# The class holding the two-point correlation functions
#
# Two-point functions have dependence on the following attributes:
#   the Momentum
#   the source time, t0
#   the source and sink operators
#   the rows of the source-sink operators
#
# For each value of the momentum, the software will average over all t0s, source/sink operators and rows.
#
class TwoPointCorrelator():
    def __init__(self, dataInfo, analysisInfo):
        self.dataInfo = dataInfo
        self.analysisInfo = analysisInfo

        # Data containers
        # "plain" means not averaged over t0,src-snk operators, rows, or momentum, i.e. there's depedence on these attributes
        self.plainData = {}    # The data that is read/loaded
        self.plainBins = {}    # The Jackknife sampling bins of the plain data
        self.plainMean = {}    # The Jackknife mean of the plain data

        self.data = {}     # The averaged data
        self.bins = {}     # The Jackknife sampling bins of the averaged data
        self.mean = {}     # The Jackknife mean of the averaged data
        
        self.Nbins = 0     # The number of Jackknife bins for the averaged data

        self.covMean = {}  # Average over all attributes but t0, needed for the Covariant Matrix

        self.binsize = self.analysisInfo['Binsize']

        self.getKey = TwoPointKeyGen

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

            for attr in ['t0','Ncfg','Nt','Nrows','Compute X-rows']:
                self.dSetAttr[mTag][attr] = dSet[attr]


            # Read source-sink operators
            intOpFile = dSet['Interpolating Operators File']
            
            self.dSetAttr[mTag]['intOpList'] = []

            with open(intOpFile) as fp:
                ops = fp.readlines()
                for op in ops:
                    self.dSetAttr[mTag]['intOpList'].append((op.split()[0],op.split()[1]))
            self.dSetAttr[mTag]['Nop'] = len(ops)

    # End __init__() -------------

    def printInfo(self):
        JSONio.dumpDictObject(self.dataInfo, '\nTwoPointCorrelator - Got the following 2pt data Info:')

        print('\nParsed the following momenta:')
        print(self.moms)

        JSONio.dumpDictObject(self.dSetAttr, '\nTwoPointCorrelator - Parsed the following Attributes:')        
    #-------------------------------

    def getData(self):
        dataSource = self.dataInfo['Data Source']
        
        if dataSource == 'ASCII':
            print('\nWill read data from ASCII files')

            for mom in self.moms:
                mTag = tags.momString(mom)
                mFTag = tags.momFile(mom)
                t0List = self.dSetAttr[mTag]['t0']
                Nrows = self.dSetAttr[mTag]['Nrows']

                # These are the dimensions of each dataset
                Ncfg = self.dSetAttr[mTag]['Ncfg']
                Nt = self.dSetAttr[mTag]['Nt']

                self.plainData[mTag] = {}
                print('Reading two-point data for momentum %s from files:'%(mTag))

                for it0,t0 in enumerate(t0List):
                    t0Tag = tags.t0(t0)
                    fileDir = ioForm.getTwoPointDirASCII(self.dataInfo['Data Main Directory'],t0Tag,mFTag)

                    for iop,opPair in enumerate(self.dSetAttr[mTag]['intOpList']):
                        srcOp,snkOp = opPair

                        for ir,row in enumerate(range(1,Nrows+1)):
                            dkey = (t0,iop,row)

                            fileName = ioForm.getTwoPointFileNameASCII(self.phaseTag,t0Tag,srcOp,snkOp,row,mFTag,self.Nvec)
                            fileRead = '%s/%s'%(fileDir,fileName)

                            print(fileRead)
                            
                            self.plainData[mTag][dkey] = np.zeros((Ncfg,Nt), dtype=np.complex128)
                            
                            with open(fileRead) as fp:
                                line = fp.readlines()
                                c = 0
                                for n in line:
                                    it   = c%Nt
                                    icfg = c//Nt
                                    self.plainData[mTag][dkey][icfg,it] = complex(np.float64(n.split()[1]),
                                                                                    np.float64(n.split()[2]))
                                    c += 1
                            
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

            Nrows = self.dSetAttr[mTag]['Nrows']
            Nt0 = len(t0List)
            Nop = self.dSetAttr[mTag]['Nop']
            Ncfg = self.dSetAttr[mTag]['Ncfg']
            Nt = self.dSetAttr[mTag]['Nt']

            Navg = Nrows + Nt0 + Nop

            # Determine the Jackknife sampling number of Bins
            self.Nbins = jackknife.Nbins(Ncfg,self.binsize)

            self.plainBins[mTag] = {}
            self.plainMean[mTag] = {}
            
            # That's the averaged data
            self.data[mTag] = np.zeros((Ncfg,Nt), dtype=np.complex128) 

            for it0,t0 in enumerate(t0List):
                for iop,opPair in enumerate(self.dSetAttr[mTag]['intOpList']):
                    for ir,row in enumerate(range(1,Nrows+1)):
                        dkey = (t0,iop,row)

                        # Jackknife sampling on the Plain data
                        self.plainBins[mTag][dkey] = np.zeros((self.Nbins,Nt), dtype=np.complex128)
                        for t in range(Nt):
                            self.plainBins[mTag][dkey][:,t] = jackknife.sampling(self.plainData[mTag][dkey][:,t].real, self.Nbins, self.binsize)

                        self.plainMean[mTag][dkey] = jackknife.mean(self.plainBins[mTag][dkey], self.Nbins, Nspl=Nt)

                        # Sum over Source-Sink operators, t0's and rows
                        self.data[mTag] += self.plainData[mTag][dkey]


            # Sum over Source-Sink operators, t0's and rows
            self.data[mTag] = self.data[mTag] / Navg

            # Jackknife sampling over the averaged data, for each momentum
            self.bins[mTag] = np.zeros((self.Nbins,Nt), dtype=np.complex128)
            for t in range(Nt):
                self.bins[mTag][:,t] = jackknife.sampling(self.data[mTag][:,t].real, self.Nbins, self.binsize)

            self.mean[mTag] = jackknife.mean(self.bins[mTag], self.Nbins, Nspl=Nt)

        print('Statistical evaluation completed')
    # End doStatistics() -------------
