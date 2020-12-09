'''
Created on Dec. 8, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Class definition that holds the ratio of three- to two-point functions
'''

import pymela.io.json_io as JSONio
import pymela.io.file_formats as ioForm
import pymela.tools.tag_creators as tags
import pymela.tools.jackknife as jackknife
import pymela.tools.gamma as gmat

import numpy as np
import h5py


# The class holding the ratio of three- to two-point functions
#
class ThreeToTwoPointCorrRatio():
    def __init__(self, c2pt, c3pt, analysisInfo):
        self.c2pt = c2pt
        self.c3pt = c3pt

        self.analysisInfo = analysisInfo

        #   The three-point function also has a significant real and imaginary part
        self.RI = ['Re','Im']

        # Ratio types
        self.ratioTypes = ['plain','sum','r-sum']

        # Data containers        
        self.bins = {}
        self.mean = {}
        for t in self.ratioTypes:
            self.bins[t] = {}
            self.mean[t] = {}
            for ri in self.RI:
                self.bins[t][ri] = {}
                self.mean[t][ri] = {}


        # Make some checks to ensure compatibility between two- and three-point functions
        if self.c2pt.Nbins != self.c3pt.Nbins:
            raise ValueError('\n Two- and three-point functions must have the same number of Jackknife bins!')

        if self.c2pt.momAvg != self.c3pt.momAvg:
            raise ValueError('\n Two- and three-point functions must have the same momenta!')

        self.momAvg = self.c2pt.momAvg
        self.Nbins = self.c2pt.Nbins

        self.dispAvg = self.c3pt.dispAvg

        self.dSetAttr3pt = self.c3pt.dSetAttr
        self.dSetAttr2pt = self.c2pt.dSetAttr
    # End __init__() -------------

    def evaluate(self):

        for mom in self.momAvg:
            mTag = tags.momString(mom)

            tsepList    = self.dSetAttr3pt[mTag]['tsep']
            tsepList_rs = self.dSetAttr3pt[mTag]['tsep'][:-1] # Need this for the reduced-summed ratio
            dispListAvg = self.dispAvg[mTag]
            gammaList   = self.dSetAttr3pt[mTag]['gamma']

            for t in self.ratioTypes:
                for ri in self.RI:
                    self.bins[t][ri][mTag] = {}
                    self.mean[t][ri][mTag] = {}

            for its, tsep in enumerate(tsepList):
                Ntins = tsep
                for z3 in dispListAvg:
                    for gamma in gammaList:
                        dkey = (tsep,z3,gamma)

                        for ri in self.RI:
                            self.bins['plain'][ri][mTag][dkey] = np.zeros((self.Nbins,Ntins),dtype = np.float128)
                            self.bins['sum'][ri][mTag][dkey]   = np.zeros(self.Nbins,dtype = np.float128)

                            for tins in range(Ntins):
                                # Plain ratio
                                self.bins['plain'][ri][mTag][dkey][:,tins] = (self.c3pt.bins[ri][mTag][dkey][:,tins] /
                                                                              self.c2pt.bins[mTag][:,tsep])

                                # Summed ratio
                                if tins>0: # Exclude source contact term
                                    self.bins['sum'][ri][mTag][dkey] += self.bins['plain'][ri][mTag][dkey][:,tins]

                            self.mean['plain'][ri][mTag][dkey] = jackknife.mean(self.bins['plain'][ri][mTag][dkey],self.Nbins, Nspl=Ntins)
                            self.mean['sum'][ri][mTag][dkey]   = jackknife.mean(self.bins['sum'][ri][mTag][dkey]  ,self.Nbins, Nspl=1)
            # End for tsep


            # Evaluate reduced-summed ratio
            for its, tsep in enumerate(tsepList_rs):
                tsepL = tsepList[its]
                tsepH = tsepList[its+1]
                for z3 in dispListAvg:
                    for gamma in gammaList:
                        dkey = (tsep,z3,gamma)
                        dkeyL = (tsepL,z3,gamma)
                        dkeyH = (tsepH,z3,gamma)
                        for ri in self.RI:
                            self.bins['r-sum'][ri][mTag][dkey] = ((self.bins['sum'][ri][mTag][dkeyH] - self.bins['sum'][ri][mTag][dkeyL]) / 
                                                                    (tsepH - tsepL))
                            self.mean['r-sum'][ri][mTag][dkey] = jackknife.mean(self.bins['r-sum'][ri][mTag][dkey],self.Nbins, Nspl=1)                                

            print('Ratio evaluation for %s completed'%(mTag))
        # End for momentum
    # End __evaluate__() -------------


    # def writeHDF5(self):
    #     h5_file = h5py.File(self.dataInfo['HDF5 Output File'],'w')

    #     # Write the Pz- and z3-averaged data
    #     for mom in self.momAvg:
    #         mTag = tags.momString(mom)
    #         mh5Tag = tags.momH5(mom)

    #         mTagD = tags.momString(mom)
    #         if mTag not in self.dSetAttr.keys():
    #             mTagD = tags.momString([mom[0],mom[1],-mom[2]])


    #         tsepList = self.dSetAttr[mTagD]['tsep']
    #         dispListAvg = self.dispAvg[mTag]

    #         for z3 in dispListAvg:
    #             dispTag = tags.disp(z3)
    #             for tsep in tsepList:
    #                 tsepTag = tags.tsep(tsep)
    #                 for gamma in self.dSetAttr[mTagD]['gamma']:
    #                     insTag = tags.insertion(gamma)
    #                     dkeyAvg = (tsep,z3,gamma)

    #                     # Write the averaged data
    #                     for ri in self.RI:
    #                         avg_group = 'avg/%s/%s/%s/%s/%s'%(mh5Tag,tsepTag,dispTag,insTag,ri)
    #                         dset_name_data = avg_group + '/data'
    #                         dset_name_bins = avg_group + '/bins'
    #                         dset_name_mean = avg_group + '/mean'
    #                         h5_file.create_dataset(dset_name_data, data = self.data[ri][mTag][dkeyAvg])
    #                         h5_file.create_dataset(dset_name_bins, data = self.bins[ri][mTag][dkeyAvg])
    #                         h5_file.create_dataset(dset_name_mean, data = self.mean[ri][mTag][dkeyAvg],dtype='f')
    #     #--------------------------------------

    #     for mom in self.moms:
    #         mTag = tags.momString(mom)
    #         mh5Tag = tags.momH5(mom)
    #         t0List = self.dSetAttr[mTag]['t0']
    #         tsepList = self.dSetAttr[mTag]['tsep']
    #         dispList = self.dSetAttr[mTag]['disp']
    #         Nrows = self.dSetAttr[mTag]['Nrows']

    #         for z3 in dispList:
    #             dispTag = tags.disp(z3)
    #             for tsep in tsepList:
    #                 tsepTag = tags.tsep(tsep)
    #                 for t0 in t0List:
    #                     t0Tag = tags.t0(t0)
    #                     for iop,opPair in enumerate(self.dSetAttr[mTag]['intOpList']):
    #                         opTag = tags.src_snk(opPair)
    #                         for row in range(1,Nrows+1):
    #                             rowTag = tags.row(row)
    #                             for gamma in self.dSetAttr[mTag]['gamma']:
    #                                 insTag = tags.insertion(gamma)
    #                                 dkey = (tsep,t0,z3,iop,row,gamma)

    #                                 # Write the plain data
    #                                 for ri in self.RI:
    #                                     plain_group = 'plain/%s/%s/%s/%s/%s/%s/%s/%s'%(mh5Tag,dispTag,tsepTag,t0Tag,opTag,rowTag,insTag,ri)
    #                                     dset_name_plainData = plain_group + '/data'
    #                                     dset_name_plainBins = plain_group + '/bins'
    #                                     dset_name_plainMean = plain_group + '/mean'

    #                                     h5_file.create_dataset(dset_name_plainData, data = self.plainData[ri][mTag][dkey])
    #                                     h5_file.create_dataset(dset_name_plainBins, data = self.plainBins[ri][mTag][dkey])
    #                                     h5_file.create_dataset(dset_name_plainMean, data = self.plainMean[ri][mTag][dkey],dtype='f')                                
    #     #--------------------------------------

    #     h5_file.close()
    #     print('Three-point function data written in HDF5.')




