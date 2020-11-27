'''
Created on Nov.20, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Class definition that holds data and manipulates Two-point correlation functions
'''

import io.json_io as JSONio
import io.file_formats as ioForm
import tools.tag_creators as tags

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
        
        self.bins = {}     # The Jackknife sampling bins of the averaged data
        self.mean = {}     # The Jackknife mean of the averaged data
        self.stdMean = {}  # The standard mean of the averaged data
        
        self.covMean = {}  # Average over all attributes but t0, needed for the Covariant Matrix
        
        self.getKey = TwoPointKeyGen


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
            srcOpFile = dSet['Source Operators File']
            snkOpFile = dSet['Sink Operators File']
            
            self.dSetAttr[mTag]['srcOpList'] = []
            self.dSetAttr[mTag]['snkOpList'] = []

            with open(srcOpFile) as fp:
                self.dSetAttr[mTag]['srcOpList'].append(fp.readlines()[0].split()[0])                    

            with open(snkOpFile) as fp:
                self.dSetAttr[mTag]['snkOpList'].append(fp.readlines()[0].split()[0])
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
                t0List = self.dSetAttr[mTag]['t0']
                Nrows = self.dSetAttr[mTag]['Nrows']
                Ncfg = self.dSetAttr[mTag]['Ncfg']
                Nt0 = len(t0List)
                Nt = self.dSetAttr[mTag]['Nt']

                self.plainData[mTag] = {}
                print('Reading two-point data for momentum %s from files:'%(mTag))
                for isrc,src in enumerate(self.dSetAttr[mTag]['srcOpList']):
                    for isnk,snk in enumerate(self.dSetAttr[mTag]['snkOpList']):
                        for it0,t0 in enumerate(t0List):
                            t0Tag = tags.t0(t0)
                            mFTag = tags.momFile(mom)
                            fileDir = ioForm.getTwoPointDirASCII(self.dataInfo['Data Main Directory'],t0Tag,mFTag)

                            for ir,row in enumerate(range(1,Nrows+1)):
                                dkey = (isrc,isnk,t0,row)

                                fileName = ioForm.getTwoPointFileNameASCII(self.phaseTag,t0Tag,src,snk,row,mFTag,self.Nvec)
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
                                # Reading file
                print('Reading two-point data for momentum %s completed.'%(mTag))
        else:
            raise ValueError('\nUnsupported "Data Source" = %s ' % (dataSource))
    # End getData() -------------
