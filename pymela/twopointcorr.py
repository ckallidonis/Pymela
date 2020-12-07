'''
Created on Nov.20, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Class definition that holds data and manipulates Two-point correlation functions
'''

import pymela.io.json_io as JSONio
import pymela.io.file_formats as ioForm
import pymela.tools.tag_creators as tags
import pymela.tools.jackknife as jackknife

import numpy as np
import h5py


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

        self.avgData = {}     # The averaged data
        self.avgBins = {}     # The Jackknife sampling bins of the averaged data
        self.avgMean = {}     # The Jackknife mean of the averaged data
        
        self.momData = {}     # The momentum-averaged data
        self.momBins = {}     # The Jackknife sampling bins of the momentum-averaged data
        self.momMean = {}     # The Jackknife mean of the momentum-averaged data

        self.Nbins = 0     # The number of Jackknife bins (same for plain and averaged data)

        self.covMean = {}  # Average over all attributes but t0, needed for the Covariant Matrix

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
            if momVec[0] != 0 and momVec[1] != 0:
                raise ValueError('\n Currently support non-zero momentum only in the z-direction!')                
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

        # Get the momenta that will be averaged over        
        self.momAvg = [[0,0,zm] for zm in list(dict.fromkeys(np.abs([z for x,y,z in self.moms])))]
        self.momAvg.sort()
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

                for t0 in t0List:
                    t0Tag = tags.t0(t0)
                    fileDir = ioForm.getTwoPointDirASCII(self.dataInfo['Data Main Directory'],t0Tag,mFTag)

                    for iop,opPair in enumerate(self.dSetAttr[mTag]['intOpList']):
                        srcOp,snkOp = opPair

                        for row in range(1,Nrows+1):
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

            Navg = Nrows * Nt0 * Nop

            # Determine the Jackknife sampling number of Bins
            self.Nbins = jackknife.Nbins(Ncfg,self.binsize)

            # The plain data Bins and Mean
            self.plainBins[mTag] = {}
            self.plainMean[mTag] = {}
            
            # That's the averaged data
            self.avgData[mTag] = np.zeros((Ncfg,Nt), dtype=np.complex128) 

            # The mean required for the covariant matrix
            self.covMean[mTag] = {}

            for t0 in t0List:
                covSum = np.zeros((Ncfg,Nt), dtype=np.float64)
                for iop,opPair in enumerate(self.dSetAttr[mTag]['intOpList']):
                    for row in range(1,Nrows+1):
                        dkey = (t0,iop,row)

                        # Jackknife sampling on the Plain data
                        self.plainBins[mTag][dkey] = np.zeros((self.Nbins,Nt), dtype=np.float64)
                        for t in range(Nt):
                            self.plainBins[mTag][dkey][:,t] = jackknife.sampling(self.plainData[mTag][dkey][:,t].real, self.Nbins, self.binsize)

                        self.plainMean[mTag][dkey] = jackknife.mean(self.plainBins[mTag][dkey], self.Nbins, Nspl=Nt)

                        # Sum over Source-Sink operators, t0's and rows
                        self.avgData[mTag] += self.plainData[mTag][dkey]

                        # Sum over Source-Sink operators and rows
                        covSum += self.plainData[mTag][dkey].real

                # Standard Mean and Error over source-sink operators and rows, for each t0 (for covariant matrix)
                covAvg = covSum/(Nop*Nrows) # Still a (Ncfg * Nt) array
                self.covMean[mTag][t0] = (np.mean(covAvg,axis=0),
                                          np.std(covAvg,axis=0)/np.sqrt(Ncfg)) 

            # Sum over Source-Sink operators, t0's and rows
            self.avgData[mTag] = self.avgData[mTag] / Navg

            # Jackknife sampling over the averaged data, for each momentum
            self.avgBins[mTag] = np.zeros((self.Nbins,Nt), dtype=np.float64)
            for t in range(Nt):
                self.avgBins[mTag][:,t] = jackknife.sampling(self.avgData[mTag][:,t].real, self.Nbins, self.binsize)

            self.avgMean[mTag] = jackknife.mean(self.avgBins[mTag], self.Nbins, Nspl=Nt)
        # End for momentum -------------

        # Perform averaging over momentum
        for mom in self.momAvg:
            momNeg = [mom[0],mom[1],-mom[2]]

            mTag = tags.momString(mom)
            mTagPos = mTag
            mTagNeg = tags.momString(momNeg)

            # If the positive (averaged) value does not exist in the Attribute Dictionary, take the negative
            mTagD = tags.momString(mom)
            if mTag not in self.dSetAttr.keys():
                mTagD = tags.momString(momNeg)

            Ncfg = self.dSetAttr[mTagD]['Ncfg']
            Nt = self.dSetAttr[mTagD]['Nt']

            self.momData[mTag] = np.zeros((Ncfg,Nt), dtype=np.complex128)
            self.momBins[mTag] = np.zeros((self.Nbins,Nt), dtype=np.float64)
            if mom in self.moms and momNeg in self.moms:
                self.momData[mTag] = 0.5 * (self.avgData[mTagPos] + self.avgData[mTagNeg])
                self.momBins[mTag] = 0.5 * (self.avgBins[mTagPos] + self.avgBins[mTagNeg])
            elif mom in self.moms and momNeg not in self.moms:
                self.momData[mTag] = self.avgData[mTagPos] 
                self.momBins[mTag] = self.avgBins[mTagPos]
            elif mom not in self.moms and momNeg in self.moms:
                self.momData[mTag] = self.avgData[mTagNeg] 
                self.momBins[mTag] = self.avgBins[mTagNeg]

            self.momMean[mTag] = jackknife.mean(self.momBins[mTag], self.Nbins, Nspl=Nt)

        print('Statistical evaluation completed')
    # End doStatistics() -------------

    def writeHDF5(self):
        h5_file = h5py.File(self.dataInfo['HDF5 Output File'],'w')

        for mom in self.moms:
            mTag = tags.momString(mom)
            mh5Tag = tags.momH5(mom)
            t0List = self.dSetAttr[mTag]['t0']
            Nrows = self.dSetAttr[mTag]['Nrows']

            # Write the averaged data
            avg_group = 'avg/%s'%(mh5Tag)
            dset_name_data = avg_group + '/data'
            dset_name_bins = avg_group + '/bins'
            dset_name_mean = avg_group + '/mean'

            h5_file.create_dataset(dset_name_data, data = self.avgData[mTag])
            h5_file.create_dataset(dset_name_bins, data = self.avgBins[mTag])
            h5_file.create_dataset(dset_name_mean, data = self.avgMean[mTag],dtype='f')

            for t0 in t0List:
                t0Tag = tags.t0(t0)

                # Write cov. matrix mean
                cov_group = 'cov/%s/%s'%(mh5Tag,t0Tag)
                dset_name_covMean = cov_group + '/mean'
                h5_file.create_dataset(dset_name_covMean, data = self.covMean[mTag][t0],dtype='f')

                for iop,opPair in enumerate(self.dSetAttr[mTag]['intOpList']):
                    opTag = tags.src_snk(opPair)
                    for row in range(1,Nrows+1):
                        rowTag = tags.row(row)

                        dkey = (t0,iop,row)

                        # Write the plain data
                        plain_group = 'plain/%s/%s/%s/%s'%(mh5Tag,t0Tag,opTag,rowTag)
                        dset_name_plainData = plain_group + '/data'
                        dset_name_plainBins = plain_group + '/bins'
                        dset_name_plainMean = plain_group + '/mean'

                        h5_file.create_dataset(dset_name_plainData, data = self.plainData[mTag][dkey])
                        h5_file.create_dataset(dset_name_plainBins, data = self.plainBins[mTag][dkey])
                        h5_file.create_dataset(dset_name_plainMean, data = self.plainMean[mTag][dkey],dtype='f')                                


        # Write the momentum-averaged data
        for mom in self.momAvg:
            mTag = tags.momString(mom)
            mh5Tag = tags.momH5(mom)

            momAvg_group = 'momAvg/%s'%(mh5Tag)
            dset_name_momData = momAvg_group + '/data'
            dset_name_momBins = momAvg_group + '/bins'
            dset_name_momMean = momAvg_group + '/mean'

            h5_file.create_dataset(dset_name_momData, data = self.momData[mTag])
            h5_file.create_dataset(dset_name_momBins, data = self.momBins[mTag])
            h5_file.create_dataset(dset_name_momMean, data = self.momMean[mTag],dtype='f')
        #--------------------------------

        h5_file.close()
        print('Two-point function data written in HDF5.')
    # End writeHDF5() -------------
