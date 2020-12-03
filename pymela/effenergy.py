'''
Created on Dec.3, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Class definition that holds data and manipulates Effective Energy, computed from the two-point functions
'''

import pymela.io.json_io as JSONio
import pymela.io.file_formats as ioForm
import pymela.tools.tag_creators as tags
import pymela.tools.jackknife as jackknife

import numpy as np
import h5py

import warnings
warnings.filterwarnings("error")

# The class holding the Effective Energy
#
# This class takes the two-point function object as input
class EffectiveEnergy():
    def __init__(self, c2pt, dataInfo):
        self.c2pt = c2pt

        self.dataInfo = dataInfo

        # Data containers
        # "plain" means not averaged over t0,src-snk operators, rows, or momentum, i.e. there's depedence on these attributes
        self.plainBins = {}    # The Jackknife sampling bins of the plain data
        self.plainMean = {}    # The Jackknife mean of the plain data

        self.avgBins = {}     # The Jackknife sampling bins of the averaged data
        self.avgMean = {}     # The Jackknife mean of the averaged data
        
        self.momBins = {}     # The Jackknife sampling bins of the momentum-averaged data
        self.momMean = {}     # The Jackknife mean of the momentum-averaged data


        # Get attributes from the two-point correlator
        self.Nbins = self.c2pt.Nbins  # The number of Jackknife bins (same for plain and averaged data)
        self.binsize = self.c2pt.analysisInfo['Binsize']

        self.moms = self.c2pt.moms
        self.dSetAttr = self.c2pt.dSetAttr
        self.momAvg = self.c2pt.momAvg
    # End __init__() -------------


    def compute(self):

        def logRatio(c2ptBins):
            Nt = self.dSetAttr[mTag]['Nt']
            binsArr = np.zeros((self.Nbins,Nt),dtype=np.float64)
            # Need to check element by element to avoid negative log warnings
            for b in range(self.Nbins):
                for t in range(Nt):
                    try:
                        binsArr[b,t] = np.log(c2ptBins[b,t] / c2ptBins[b,(t+1)%Nt])
                    except RuntimeWarning:
                        binsArr[b,t] = None
            meanArr = jackknife.mean(binsArr, self.Nbins, Nspl=Nt)

            return binsArr,meanArr
        #---------------------

        for mom in self.moms:
            mTag = tags.momString(mom)
            mFTag = tags.momFile(mom)
            t0List = self.dSetAttr[mTag]['t0']
            Nrows = self.dSetAttr[mTag]['Nrows']
            Nt = self.dSetAttr[mTag]['Nt']

            # Effective Energy for averaged data
            self.avgBins[mTag], self.avgMean[mTag] = logRatio(self.c2pt.avgBins[mTag])

            self.plainBins[mTag] = {}
            self.plainMean[mTag] = {}
            for it0,t0 in enumerate(t0List):
                for iop,opPair in enumerate(self.dSetAttr[mTag]['intOpList']):
                    srcOp,snkOp = opPair
                    for ir,row in enumerate(range(1,Nrows+1)):
                        dkey = (t0,iop,row)

                        # Effective Energy for plain data
                        self.plainBins[mTag][dkey], self.plainMean[mTag][dkey] = logRatio(self.c2pt.plainBins[mTag][dkey])

            # Effective Energy for momentum-averaged data
            if mTag in self.momAvg.keys():
                self.momBins[mTag], self.momMean[mTag] = logRatio(self.c2pt.momBins[mTag])

        print('Effective Energy computed.')
    # End compute() -------------

    def writeHDF5(self):
        h5_file = h5py.File(self.dataInfo['HDF5 Output File'],'w')

        for mom in self.moms:
            mTag = tags.momString(mom)
            mh5Tag = tags.momH5(mom)
            t0List = self.dSetAttr[mTag]['t0']
            Nrows = self.dSetAttr[mTag]['Nrows']

            # Write the averaged data
            avg_group = 'avg/%s'%(mh5Tag)
            dset_name_bins = avg_group + '/bins'
            dset_name_mean = avg_group + '/mean'

            dset_bins = h5_file.create_dataset(dset_name_bins, data = self.avgBins[mTag])
            dset_mean = h5_file.create_dataset(dset_name_mean, data = self.avgMean[mTag],dtype='f')

            for it0,t0 in enumerate(t0List):
                t0Tag = tags.t0(t0)

                for iop,opPair in enumerate(self.dSetAttr[mTag]['intOpList']):
                    opTag = tags.src_snk(opPair)
                    for ir,row in enumerate(range(1,Nrows+1)):
                        rowTag = tags.row(row)

                        dkey = (t0,iop,row)

                        # Write the plain data
                        plain_group = 'plain/%s/%s/%s/%s'%(mh5Tag,t0Tag,opTag,rowTag)
                        dset_name_plainBins = plain_group + '/bins'
                        dset_name_plainMean = plain_group + '/mean'

                        dset_plainBins = h5_file.create_dataset(dset_name_plainBins, data = self.plainBins[mTag][dkey])
                        dset_plainMean = h5_file.create_dataset(dset_name_plainMean, data = self.plainMean[mTag][dkey],dtype='f')                                


        # Write the momentum-averaged data
        for mKey in self.momAvg.keys():
            mh5Tag = tags.momH5(tags.momVec(mKey))

            momAvg_group = 'momAvg/%s'%(mh5Tag)
            dset_name_momBins = momAvg_group + '/bins'
            dset_name_momMean = momAvg_group + '/mean'

            dset_momBins = h5_file.create_dataset(dset_name_momBins, data = self.momBins[mKey])
            dset_momMean = h5_file.create_dataset(dset_name_momMean, data = self.momMean[mKey],dtype='f')
        #--------------------------------

        h5_file.close()
        print('Effective Energy data written in HDF5.')
    # End writeHDF5() -------------
