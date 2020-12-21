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

import numpy as np
import h5py


# The class holding the ratio of three- to two-point functions
#
class ThreeToTwoPointCorrRatio():
    def __init__(self, c2pt, c3pt, dataInfo, analysisInfo):
        self.c2pt = c2pt
        self.c3pt = c3pt

        self.dataInfo = dataInfo
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

        self.gammaList   = self.c3pt.gammaList

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

            for t in self.ratioTypes:
                for ri in self.RI:
                    self.bins[t][ri][mTag] = {}
                    self.mean[t][ri][mTag] = {}

            for its, tsep in enumerate(tsepList):
                Ntins = tsep
                for z3 in dispListAvg:
                    for gamma in self.gammaList:
                        dkey = (tsep,z3,gamma)

                        for ri in self.RI:
                            self.bins['plain'][ri][mTag][dkey] = np.zeros((self.Nbins,Ntins),dtype = np.float64)
                            self.bins['sum'][ri][mTag][dkey]   = np.zeros(self.Nbins,dtype = np.float64)

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
                    for gamma in self.gammaList:
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


    def writeHDF5(self):
        h5_file = h5py.File(self.dataInfo['HDF5 Output File'],'w')


        for mom in self.momAvg:
            mTag = tags.momString(mom)
            mh5Tag = tags.momH5(mom)            

            tsepList    = self.dSetAttr3pt[mTag]['tsep']
            tsepList_rs = self.dSetAttr3pt[mTag]['tsep'][:-1]
            dispListAvg = self.dispAvg[mTag]
            Ntsep    = len(tsepList)
            Ntsep_rs = len(tsepList_rs)

            for z3 in dispListAvg:
                dispTag = tags.disp(z3)
                for gamma in self.gammaList:
                    insTag = tags.insertion(gamma)
                    for ri in self.RI:
                        
                        sumRatioH5 = (np.zeros(Ntsep),np.zeros(Ntsep),np.zeros(Ntsep))
                        for its,tsep in enumerate(tsepList):
                            dkey = (tsep,z3,gamma)
                            tsepTag = tags.tsep(tsep)

                            # Write the plain ratio bins and mean
                            rType = 'plain'
                            group = '%s/%s/%s/%s/%s/%s'%(rType,mh5Tag,tsepTag,dispTag,insTag,ri)
                            dset_name_bins = 'bins/' + group 
                            dset_name_mean = 'mean/' + group 

                            h5_file.create_dataset(dset_name_bins, data = self.bins[rType][ri][mTag][dkey])
                            h5_file.create_dataset(dset_name_mean, data = self.mean[rType][ri][mTag][dkey],dtype='f')
                            #---------------------------------------------------------------

                            # Write the summed ratio bins
                            rType = 'sum'
                            group = '%s/%s/%s/%s/%s/%s'%(rType,mh5Tag,tsepTag,dispTag,insTag,ri)
                            dset_name_bins = 'bins/' + group
                            h5_file.create_dataset(dset_name_bins, data = self.bins[rType][ri][mTag][dkey])

                            # Convert the summed ratio mean into arrays that depend on tsep
                            sumRatioH5[0][its] = tsep # tsep (x)
                            sumRatioH5[1][its] = self.mean[rType][ri][mTag][dkey][0] # ratio mean  (y)
                            sumRatioH5[2][its] = self.mean[rType][ri][mTag][dkey][1] # ratio error (y-error)
                        # End for tsep

                        # Write the summed ratio means
                        rType = 'sum'
                        group = '%s/%s/%s/%s/%s'%(rType,mh5Tag,dispTag,insTag,ri)
                        dset_name_mean = 'mean/' + group
                        h5_file.create_dataset(dset_name_mean, data = sumRatioH5, dtype='f')
                        #-----------------------------


                        # Reduced-summed ratio
                        rSumRatioH5 = (np.zeros(Ntsep_rs),np.zeros(Ntsep_rs),np.zeros(Ntsep_rs))
                        for its,tsep in enumerate(tsepList_rs):
                            dkey = (tsep,z3,gamma)
                            tsepTag = tags.tsep(tsep)

                            # Write the summed ratio bins
                            rType = 'r-sum'
                            group = '%s/%s/%s/%s/%s/%s'%(rType,mh5Tag,tsepTag,dispTag,insTag,ri)
                            dset_name_bins = 'bins/' + group
                            h5_file.create_dataset(dset_name_bins, data = self.bins[rType][ri][mTag][dkey])

                            # Convert the reduced-summed ratio mean into arrays that depend on tsep
                            rSumRatioH5[0][its] = tsep # tsep (x)
                            rSumRatioH5[1][its] = self.mean[rType][ri][mTag][dkey][0] # ratio mean  (y)
                            rSumRatioH5[2][its] = self.mean[rType][ri][mTag][dkey][1] # ratio error (y-error)
                        # End for tsep

                        # Write the reduced-summed ratio means
                        rType = 'r-sum'
                        group = '%s/%s/%s/%s/%s'%(rType,mh5Tag,dispTag,insTag,ri)
                        dset_name_mean = 'mean/' + group
                        h5_file.create_dataset(dset_name_mean, data = rSumRatioH5, dtype='f')
                       #-----------------------------
        # End for momentum

        h5_file.close()
        print('Ratio data written in HDF5.')




