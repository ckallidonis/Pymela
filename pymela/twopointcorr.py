'''
Created on Nov.20, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Class definition that holds data and manipulates Two-point correlation functions
'''

import io.json_io as JSONio
import io.file_formats as ioForm
import tools.tag_creators as tags

def TwoPointKeyGen(mom,t0=None,isrcOp=None,isnkOp=None,iRow=None):
    mTag = tags.momentum(mom)
    if t0==None and isrcOp==None and isnkOp==None and iRow==None:
        return (mTag)
    else:
        return (mTag,t0,isrcOp,isnkOp,iRow)
#-------------------------------

    
# The class holding the two-point correlation functions
# Two-point functions have dependence on the following attributes:
#   the source time, t0
#   the source and sink operators
#   the rows of the source-sink operators
#   the Momentum
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


        # Fill in momentum and t0 lists
        self.moms  = []
        self.t0 = {}
        for dSet in self.dataInfo['Datasets']:
            try:
                momVec = dSet['mom']
            except:
                raise ValueError('Each dataset must have a "z-Mom" key, whose value must be an integer.')
            self.moms.append(momVec)

            mTag = tags.momString(momVec)
            try:
                tList = dSet['t0']
            except:
                raise ValueError('Each dataset must have a "t0" key, whose value must be a list of integers.')
            self.t0[mTag] = tList


        # Read source-sink operators
        self.snkOpFile = self.dataInfo['Sink Operator File']
        self.nRows = self.dataInfo['Operator Nrows']
        if self.dataInfo['Compute X-rows']:
            raise ValueError('Does not support doing cross-rows in two-point function for now!')
        
        self.srcOp = {}
        self.snkOp = {}
        for mom in self.moms:
            momStr = tags.momString(mom)
            
            self.srcOp[momStr] = []
            try:
                srcOpFile = self.dataInfo['Source Operator File'][momStr]                
            except:
                raise ValueError('Expected key %s in "Source Operator File" object'%(momStr))
            with open(srcOpFile) as fp:
                self.srcOp[momStr].append(fp.readlines()[0].split()[0])
                
            self.snkOp[momStr] = []
            try:
                snkOpFile = self.dataInfo['Sink Operator File'][momStr]                
            except:
                raise ValueError('Expected key %s in "Sink Operator File" object'%(momStr))
            with open(snkOpFile) as fp:
                self.snkOp[momStr].append(fp.readlines()[0].split()[0])                
    # End __init__() -------------

    def printInfo(self):
        JSONio.dumpDictObject(self.dataInfo, '\nTwoPointCorrelator - Got the following 2pt data Info:')

        print('\nParsed the following momenta:')
        print(self.moms)
        print('\nParsed the following t-sources:')
        print(self.t0)
        print('\nParsed the following Source Operators:')
        print(self.srcOp)
        print('\nParsed the following Sink Operators:')
        print(self.snkOp)
    #-------------------------------

    def getData(self):
        dataSource = self.dataInfo['Data Source']
        
        if dataSource == 'ASCII':
            print('\nWill read data from ASCII')

            for mom in self.moms:
                mTag = tags.momString(mom)
                Nt0 = len(self.t0[mTag])
                for it0,t0 in enumerate(self.t0[mTag]):
                    fileDir = ioForm.getTwoPointDirASCII(self.dataInfo['Data Main Directory'], tags.t0(t0), tags.momFile(mom))

#                    for isrc,src in enumerate(self.srcOpList):
#                        for isnk,snk in enumerate(self.snkOpList):
                    
            
        else:
            raise ValueError('\nUnsupported "Data Source" = %s ' % (dataSource))
    # End getData() -------------
