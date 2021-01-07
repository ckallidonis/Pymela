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

        # The list of insertion operators we are considering
        self.gammaList = self.dataInfo['Insertion Operators']

        self.moms  = []
        self.dispAvg = {}
        self.dSetAttr = {}
        self.dSetList = self.dataInfo['Datasets']
        for dSet in self.dSetList:
            momList = dSet['Mom List']

            if dSet['Compute X-rows']:
                raise ValueError('Does not support doing cross-rows in two-point function for now!')

            for momVec in momList:
                if momVec[0] != 0 and momVec[1] != 0:
                    raise ValueError('\n Currently support non-zero momentum only in the z-direction!')
                mTag = tags.momString(momVec) # Dataset Attributes are listed for each momentum
                self.dSetAttr[mTag] = {}

                self.moms.append(momVec)

                for attr in ['t0','Ncfg','tsep','disp','Nrows','Compute X-rows','Phase Info']:
                    self.dSetAttr[mTag][attr] = dSet[attr]

                # Determine the values of z3 that we will average over
                self.dispAvg[mTag] = list(dict.fromkeys(np.abs(self.dSetAttr[mTag]['disp'])))
                self.dispAvg[mTag].sort()

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

        # Make sure to include entries for the averaged momentum in the dataset attributes
        for mom in self.momAvg:
            mTag = tags.momString(mom)
            if mTag not in self.dSetAttr.keys():
                mTagD = tags.momString([mom[0],mom[1],-mom[2]])
                self.dSetAttr[mTag] = self.dSetAttr[mTagD]
                print('Momentum %s not in original dataset attributes. Adding from momentum %s'%(mTag,mTagD))
    # End __init__() -------------

    def printInfo(self):
        JSONio.dumpDictObject(self.dataInfo, '\nThreePointCorrelator - Got the following 3pt data Info:')

        print('\nParsed the following momenta:')
        print(self.moms)

        JSONio.dumpDictObject(self.dSetAttr, '\nThreePointCorrelator - Parsed the following Attributes:')        

        print('\n Will average over the momenta:', self.momAvg)
        print('\n Will average over the displacement values for each momentum:', self.dispAvg)
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
                phaseInfo = self.dSetAttr[mTag]['Phase Info']

                # Determine phase tag based on momentum sign
                if list(phaseInfo.keys())[0] == 'unphased':
                    phFile = phaseInfo['unphased']
                    phDir = 'unphased'
                elif list(phaseInfo.keys())[0] == 'phased':
                    phFile = phaseInfo['phased']['Plus'] if mom[2] >= 0 else phaseInfo['phased']['Minus']
                    phDir = 'phased/' + phFile
                else:
                    raise ValueError('Supported Phase Tag keys are ["unphased","phased"]')

                # Determine the Jackknife sampling number of Bins
                self.Nbins = jackknife.Nbins(Ncfg,self.binsize)

                for ri in self.RI:
                    self.plainData[ri][mTag] = {}

                for tsep in tsepList:
                    Nt = tsep
                    tsepTag = tags.tsep(tsep)
                    for t0 in t0List:
                        t0Tag = tags.t0(t0)
                        fileDir = ioForm.getThreePointDirASCII(self.dataInfo['Data Main Directory'],phDir,t0Tag,tsepTag,mFTag)
                        print('Reading three-point data for momentum %s, tsep = %d, t0 = %d'%(mTag,tsep,t0))
                        for z3 in dispList:
                            dispTag = tags.disp(z3)
                            for iop,opPair in enumerate(self.dSetAttr[mTag]['intOpList']):
                                srcOp,snkOp = opPair
                                for row in range(1,Nrows+1):

                                    for gamma in self.gammaList:
                                        dkey = (tsep,t0,z3,iop,row,gamma)

                                        # Determine gamma matrix name and row
                                        insOp,insRow = gmat.insertionMap(gamma)

                                        fileName = ioForm.getThreePointFileNameASCII(phFile,t0Tag,tsepTag,
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
                                                rawData[icfg,it] = complex(np.float128(n.split()[1]),
                                                                           np.float128(n.split()[2]))
                                                c += 1
                                        # Done reading file

                                        self.plainData['Re'][mTag][dkey] = rawData.real
                                        self.plainData['Im'][mTag][dkey] = rawData.imag

                print('Reading three-point data for momentum %s completed.\n'%(mTag))

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
                    for gamma in self.gammaList:
                        dkeyAvg = (tsep,z3,gamma)
                        for ri in self.RI:
                            self.avgData[ri][mTag][dkeyAvg] = np.zeros((Ncfg,Nt),dtype=np.float128)

                        # We are averaging for the following attributes
                        for t0 in t0List:
                            for iop,opPair in enumerate(self.dSetAttr[mTag]['intOpList']):
                                for row in range(1,Nrows+1):
                                    dkey = (tsep,t0,z3,iop,row,gamma)

                                    for ri in self.RI:
                                        # Jackknife sampling on the Plain data
                                        self.plainBins[ri][mTag][dkey] = np.zeros((self.Nbins,Nt), dtype=np.float128)
                                        for t in range(Nt):
                                            self.plainBins[ri][mTag][dkey][:,t] = jackknife.sampling(self.plainData[ri][mTag][dkey][:,t], self.Nbins, self.binsize)
                                        self.plainMean[ri][mTag][dkey] = jackknife.mean(self.plainBins[ri][mTag][dkey], self.Nbins, Nspl=Nt)

                                        # Average over Source-Sink operators, t0's and rows
                                        self.avgData[ri][mTag][dkeyAvg] += self.plainData[ri][mTag][dkey]

                        # Average over Source-Sink operators, t0's and rows
                        for ri in self.RI:
                            self.avgData[ri][mTag][dkeyAvg] = self.avgData[ri][mTag][dkeyAvg] / Navg

                            # Jackknife sampling over the averaged data, for each momentum, tsep, z3 and gamma
                            self.avgBins[ri][mTag][dkeyAvg] = np.zeros((self.Nbins,Nt), dtype=np.float128)
                            for t in range(Nt):
                                self.avgBins[ri][mTag][dkeyAvg][:,t] = jackknife.sampling(self.avgData[ri][mTag][dkeyAvg][:,t],
                                                                                          self.Nbins, self.binsize)
                            self.avgMean[ri][mTag][dkeyAvg] = jackknife.mean(self.avgBins[ri][mTag][dkeyAvg],
                                                                             self.Nbins, Nspl=Nt)

            print('Jackknife analysis for momentum %s completed'%(mTag))
        # End for momentum


        # Perform average over momenta and z3 values
        for mom in self.momAvg:
            mTag = tags.momString(mom)

            tsepList = self.dSetAttr[mTag]['tsep']
            dispList = self.dSetAttr[mTag]['disp']
            dispListAvg = self.dispAvg[mTag]
            Ncfg = self.dSetAttr[mTag]['Ncfg']

            for ri in self.RI:
                self.data[ri][mTag] = {}
                self.bins[ri][mTag] = {}
                self.mean[ri][mTag] = {}

            for tsep in tsepList:
                Nt = tsep
                for gamma in self.gammaList:
                    for z3 in dispListAvg: # Run over the z3>=0
                        dkey = (tsep,z3,gamma)

                        for ri in self.RI:
                            self.data[ri][mTag][dkey] = np.zeros((Ncfg,Nt), dtype=np.float128)

                        if mom == [0,0,0]:
                            if z3 == 0 or not (z3 in dispList and -z3 in dispList): # Pz=0, z3=0, OR NOT both z3 and -z3 exist
                                dkeyAvg = (tsep,-z3,gamma) if -z3 in dispList else (tsep,z3,gamma)
                                for ri in self.RI:
                                    self.data[ri][mTag][dkey] = self.avgData[ri][mTag][dkeyAvg]
                            else: # Pz=0, z3!=0
                                dkeyAvgPosZ = (tsep, z3,gamma)
                                dkeyAvgNegZ = (tsep,-z3,gamma)
                                self.data['Re'][mTag][dkey] = 0.5 * (self.avgData['Re'][mTag][dkeyAvgPosZ] + self.avgData['Re'][mTag][dkeyAvgNegZ])
                                self.data['Im'][mTag][dkey] = 0.5 * (self.avgData['Im'][mTag][dkeyAvgPosZ] - self.avgData['Im'][mTag][dkeyAvgNegZ])
                        else:
                            momNeg = [mom[0],mom[1],-mom[2]]
                            if mom in self.moms and momNeg in self.moms: # Negative momentum exists in global momentum list
                                mTagPos = mTag
                                mTagNeg = tags.momString(momNeg)
                                if z3 == 0: # Pz!=0, z3=0
                                    self.data['Re'][mTag][dkey] = 0.5 * (self.avgData['Re'][mTagPos][dkey] + self.avgData['Re'][mTagNeg][dkey])
                                    self.data['Im'][mTag][dkey] = 0.5 * (self.avgData['Im'][mTagPos][dkey] - self.avgData['Im'][mTagNeg][dkey])
                                else: # Pz!=0, z3!=0
                                    if z3 in dispList and -z3 in dispList:
                                        dkeyAvgPosZ = (tsep, z3,gamma)
                                        dkeyAvgNegZ = (tsep,-z3,gamma)

                                        self.data['Re'][mTag][dkey] = 0.25 * (self.avgData['Re'][mTagPos][dkeyAvgPosZ] +
                                                                              self.avgData['Re'][mTagPos][dkeyAvgNegZ] +
                                                                              self.avgData['Re'][mTagNeg][dkeyAvgPosZ] +
                                                                              self.avgData['Re'][mTagNeg][dkeyAvgNegZ])

                                        self.data['Im'][mTag][dkey] = 0.25 * (self.avgData['Im'][mTagPos][dkeyAvgPosZ] -
                                                                              self.avgData['Im'][mTagPos][dkeyAvgNegZ] -
                                                                              self.avgData['Im'][mTagNeg][dkeyAvgPosZ] +
                                                                              self.avgData['Im'][mTagNeg][dkeyAvgNegZ])
                                    elif z3 in dispList and -z3 not in dispList:
                                        dkeyAvg = (tsep,z3,gamma)
                                        self.data['Re'][mTag][dkey] = 0.5 * (self.avgData['Re'][mTagPos][dkeyAvg] +
                                                                             self.avgData['Re'][mTagNeg][dkeyAvg])
                                        self.data['Im'][mTag][dkey] = 0.5 * (self.avgData['Im'][mTagPos][dkeyAvg] -
                                                                             self.avgData['Im'][mTagNeg][dkeyAvg])
                                    elif -z3 in dispList and z3 not in dispList:
                                        dkeyAvg = (tsep,-z3,gamma)
                                        self.data['Re'][mTag][dkey] = 0.5 * (self.avgData['Re'][mTagPos][dkeyAvg] +
                                                                             self.avgData['Re'][mTagNeg][dkeyAvg])
                                        self.data['Im'][mTag][dkey] = 0.5 * (self.avgData['Im'][mTagNeg][dkeyAvg] -
                                                                             self.avgData['Im'][mTagPos][dkeyAvg])
                                    else:
                                        raise ValueError('\n Error: Inconsistency with z3 values!!!')
                            elif mom in self.moms and momNeg not in self.moms:
                                mTagPos = mTag
                                if z3 == 0 or not (z3 in dispList and -z3 in dispList): # Pz!=0, z3=0
                                    dkeyAvg = (tsep,-z3,gamma) if -z3 in dispList else (tsep,z3,gamma)
                                    for ri in self.RI:
                                        self.data[ri][mTag][dkey] = self.avgData[ri][mTagPos][dkeyAvg]
                                else: # Pz!=0, z3!=0
                                    dkeyAvgPosZ = (tsep, z3,gamma)
                                    dkeyAvgNegZ = (tsep,-z3,gamma)

                                    self.data['Re'][mTag][dkey] = 0.5 * (self.avgData['Re'][mTagPos][dkeyAvgPosZ] +
                                                                         self.avgData['Re'][mTagPos][dkeyAvgNegZ])

                                    self.data['Im'][mTag][dkey] = 0.5 * (self.avgData['Im'][mTagPos][dkeyAvgPosZ] -
                                                                         self.avgData['Im'][mTagPos][dkeyAvgNegZ])
                            elif momNeg in self.moms and mom not in self.moms:
                                mTagNeg = tags.momString(momNeg)
                                if z3 == 0 or not (z3 in dispList and -z3 in dispList): # Pz!=0, z3=0
                                    dkeyAvg = (tsep,-z3,gamma) if -z3 in dispList else (tsep,z3,gamma)
                                    for ri in self.RI:
                                        self.data[ri][mTag][dkey] = self.avgData[ri][mTagNeg][dkeyAvg]
                                else: # Pz!=0, z3!=0
                                    dkeyAvgPosZ = (tsep, z3,gamma)
                                    dkeyAvgNegZ = (tsep,-z3,gamma)

                                    self.data['Re'][mTag][dkey] = 0.5 * (self.avgData['Re'][mTagNeg][dkeyAvgPosZ] +
                                                                         self.avgData['Re'][mTagNeg][dkeyAvgNegZ])

                                    self.data['Im'][mTag][dkey] = 0.5 * (self.avgData['Im'][mTagNeg][dkeyAvgNegZ] -
                                                                         self.avgData['Im'][mTagNeg][dkeyAvgPosZ])
                            else:
                                raise ValueError('\n Error: Inconsistency with momenta values!!!')
                        # End if mom != 0

                        # Jackknife sampling over the fully averaged data, for each momentum, tsep, z3 and gamma
                        for ri in self.RI:
                            self.bins[ri][mTag][dkey] = np.zeros((self.Nbins,Nt), dtype=np.float128)
                            for t in range(Nt):
                                self.bins[ri][mTag][dkey][:,t] = jackknife.sampling(self.data[ri][mTag][dkey][:,t], self.Nbins, self.binsize)
                            self.mean[ri][mTag][dkey] = jackknife.mean(self.bins[ri][mTag][dkey], self.Nbins, Nspl=Nt)

            print('Averaging over z3 and momenta for momentum %s completed.'%(mTag))

    def writeHDF5(self):
        h5_file = h5py.File(self.dataInfo['HDF5 Output File'],'w')

        # Write the Pz- and z3-averaged data
        for mom in self.momAvg:
            mTag = tags.momString(mom)
            mh5Tag = tags.momH5(mom)

            tsepList = self.dSetAttr[mTag]['tsep']
            dispListAvg = self.dispAvg[mTag]

            for z3 in dispListAvg:
                dispTag = tags.disp(z3)
                for tsep in tsepList:
                    tsepTag = tags.tsep(tsep)
                    for gamma in self.gammaList:
                        insTag = tags.insertion(gamma)
                        dkeyAvg = (tsep,z3,gamma)

                        # Write the averaged data
                        for ri in self.RI:
                            avg_group = 'avg/%s/%s/%s/%s/%s'%(mh5Tag,tsepTag,dispTag,insTag,ri)
                            dset_name_data = avg_group + '/data'
                            dset_name_bins = avg_group + '/bins'
                            dset_name_mean = avg_group + '/mean'
                            h5_file.create_dataset(dset_name_data, data = self.data[ri][mTag][dkeyAvg])
                            h5_file.create_dataset(dset_name_bins, data = self.bins[ri][mTag][dkeyAvg])
                            h5_file.create_dataset(dset_name_mean, data = self.mean[ri][mTag][dkeyAvg],dtype='f')
        #--------------------------------------

        for mom in self.moms:
            mTag = tags.momString(mom)
            mh5Tag = tags.momH5(mom)
            t0List = self.dSetAttr[mTag]['t0']
            tsepList = self.dSetAttr[mTag]['tsep']
            dispList = self.dSetAttr[mTag]['disp']
            Nrows = self.dSetAttr[mTag]['Nrows']

            for z3 in dispList:
                dispTag = tags.disp(z3)
                for tsep in tsepList:
                    tsepTag = tags.tsep(tsep)
                    for t0 in t0List:
                        t0Tag = tags.t0(t0)
                        for iop,opPair in enumerate(self.dSetAttr[mTag]['intOpList']):
                            opTag = tags.src_snk(opPair)
                            for row in range(1,Nrows+1):
                                rowTag = tags.row(row)
                                for gamma in self.gammaList:
                                    insTag = tags.insertion(gamma)
                                    dkey = (tsep,t0,z3,iop,row,gamma)

                                    # Write the plain data
                                    for ri in self.RI:
                                        plain_group = 'plain/%s/%s/%s/%s/%s/%s/%s/%s'%(mh5Tag,dispTag,tsepTag,t0Tag,opTag,rowTag,insTag,ri)
                                        dset_name_plainData = plain_group + '/data'
                                        dset_name_plainBins = plain_group + '/bins'
                                        dset_name_plainMean = plain_group + '/mean'

                                        h5_file.create_dataset(dset_name_plainData, data = self.plainData[ri][mTag][dkey])
                                        h5_file.create_dataset(dset_name_plainBins, data = self.plainBins[ri][mTag][dkey])
                                        h5_file.create_dataset(dset_name_plainMean, data = self.plainMean[ri][mTag][dkey],dtype='f')                                
        #--------------------------------------

        h5_file.close()
        print('Three-point function data written in HDF5.')




