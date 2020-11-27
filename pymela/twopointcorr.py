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


        # Fill in Attributes
        self.moms  = []
        self.dSetAttr = {}
        dSetList = self.dataInfo['Datasets']
        for dSet in dSetList:
            momVec = dSet['mom']
            mTag = tags.momString(momVec) # Dataset Attributes are listed for each momentum
            self.dSetAttr[mTag] = {}

            self.moms.append(momVec)

            if dSet['Compute X-rows']:
                raise ValueError('Does not support doing cross-rows in two-point function for now!')

            for attr in ['t0','Ncfg','Nt', 'Operator Nrows', 'Compute X-rows']:
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
            print('\nWill read data from ASCII')

            # for mom in self.moms:
            #     mTag = tags.momString(mom)
            #     Nt0 = len(self.t0[mTag])
            #     for it0,t0 in enumerate(self.t0[mTag]):
            #         fileDir = ioForm.getTwoPointDirASCII(self.dataInfo['Data Main Directory'], tags.t0(t0), tags.momFile(mom))

#                    for isrc,src in enumerate(self.srcOpList):
#                        for isnk,snk in enumerate(self.snkOpList):
                    
            
        else:
            raise ValueError('\nUnsupported "Data Source" = %s ' % (dataSource))
    # End getData() -------------
