'''
Created on Dec. 15, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Class definition that computes reduced Ioffe-time distributions
'''

import pymela.io.json_io as JSONio
import pymela.io.file_formats as ioForm
import pymela.tools.tag_creators as tags
import pymela.tools.jackknife as jackknife
import pymela.tools.gamma as gmat

import numpy as np
import h5py
import scipy.optimize as scipyOpt

# The class holding the Ioffe-time Distributions
#
class ITD():
    def __init__(self, plat = None, summ = None, ITDinfo = None, fitInfo = None, ensembleInfo = None):

        if plat == None and summ == None:
            raise ValueError('All of the supported fit types are "None". Cannot define ITDs!')

        self.fitInfo = fitInfo
        self.info = ITDinfo
        self.ensInfo = ensembleInfo

        # Real-Imaginary part
        self.RI = ['Re','Im']

        # Momentum-displacement position list
        # 'c'   : p = mom    , disp = z3 (Center point)
        # 'z0'  : p = mom    , disp = 0
        # 'p0'  : p = (0,0,0), disp = z3
        # 'p0z0': p = (0,0,0), disp = 0
        self.pzPos = ['c','z0','p0','p0z0']
        self.pzPosOff = ['z0','p0','p0z0'] # No center point


        # Types of fits that we are considering for the ITDs
        self.fitLabels = self.info['Optimal Fits'].keys()

        self.fitTypes = {'Plateau': [], 'Summation':[]}

        # Make sure that the input labels are included in the fits performed earlier
        self.plat = plat
        if self.plat != None:
            for fitSeq in self.plat.fitInfo:
                fLabel = fitSeq['Label']
                if fLabel not in self.fitLabels:
                    raise ValueError('Fit Label %s not in Input Fit Labels'%(fLabel))
                self.fitTypes['Plateau'].append(fLabel)

        self.summ = summ
        if self.summ != None:
            for fitSeq in self.summ.fitInfo:
                fLabel = fitSeq['Label']
                if fLabel not in self.fitLabels:
                    raise ValueError('Fit Label %s not in Input Fit Labels'%(fLabel))
                self.fitTypes['Summation'].append(fLabel)
        #--------------------------------------

        if self.plat != None: 
            # If self.summ also != None then that's still OK, because these attributes are the same
            # by definition in plat and summ
            self.momAvg  = self.plat.momAvg
            self.dispAvg = self.plat.dispAvg
            self.Nbins   = self.plat.Nbins
            self.dSetAttr3pt = self.plat.dSetAttr3pt
        else:
            # Get these attributes from the summ fits instead, it MUST be defined otherwise ValueError is raised
            self.momAvg  = self.summ.momAvg
            self.dispAvg = self.summ.dispAvg
            self.Nbins   = self.summ.Nbins
            self.dSetAttr3pt = self.summ.dSetAttr3pt


        # The ITD bins and mean
        self.bins = {}
        self.mean = {}

        # Read-in the selected fits that will be used in ITD evaluation
        self.tSelFit = {}

        for fit in self.fitLabels:
            self.bins[fit] = {}
            self.mean[fit] = {}
            self.tSelFit[fit] = {}
            for ri in self.RI:
                self.tSelFit[fit][ri] = {}

                for tOpt,momDisp in self.info['Optimal Fits'][fit][ri].items():
                    for md in momDisp:
                        for mom in md[0]:
                            mTag = tags.momString(mom)
                            for z3 in md[1]:
                                self.tSelFit[fit][ri][(mTag,z3)] = int(tOpt)
        #--------------------------


        # Determine if there are plateau ranges provided
        # Only Real-Imaginary and tsep dependence for now
        self.platRng = None
        if 'Plateau Ranges' in self.info.keys():
            self.platRng = {}
            for ri in self.RI:
                self.platRng[ri] = {}
                for t,r in self.info['Plateau Ranges'][ri].items():
                    self.platRng[ri][int(t)] = r

        print('ITD initialized')
    # End __init__() -------------


    def evaluate(self):

        # This function can be used by other fit types as well
        def getSelectedFits(fit_,ri_,mTag_c,z3_c):
            mTag_0 = tags.momString([0,0,0])
            z3_0 = 0
            t = {}
            for pz in self.pzPos:
                zT = z3_0   if 'z0' in pz else z3_c 
                mT = mTag_0 if 'p0' in pz else mTag_c
                t[pz] = self.tSelFit[fit_][ri_][(mT,zT)]
            return t
        #----------------------


        def evaluatePlateauITD(fitLabels):
            # This function is specific to the plateau fits
            def getOptimalPlatKeys(tOpt,z3_c,gamma):
                z3_0 = 0
                kp = {}
                for pz in self.pzPos:
                    zT = z3_0 if 'z0' in pz else z3_c
                    kp[pz] = (tOpt[pz],zT,gamma)
                return kp
            #-------------------

            # This function is specific to the plateau fits
            def getOptimalPlatFits(fit_,ri_,mTag_c,dkey_):
                mTag_0 = tags.momString([0,0,0])
                f = {}                
                for pz in self.pzPos:
                    mT = mTag_0 if 'p0' in pz else mTag_c
                    if self.plat.optimalFit[fit_][ri_][mT][dkey_[pz]] != -1:
                        # Valid optimal fit
                        f[pz] = self.plat.optimalFit[fit_][ri_][mT][dkey_[pz]]
                    else:
                        # Get value from provided plateau ranges
                        # These are the same for all momenta and z
                        tS = dkey_[pz][0] # Get the tsep from the key
                        f[pz] = self.platRng[ri_][tS]
                return f
            #-------------------

            # Zero momentum
            mTag_0 = tags.momString([0,0,0])

            # Temporary variable
            platBins = {}
            for ri in self.RI:
                platBins[ri] = {}

            for fit in fitLabels: # These are just the plateau fit labels!
                for mom in self.momAvg:
                    mTag = tags.momString(mom)
                    dispListAvg = self.dispAvg[mTag]
                    gammaList   = self.dSetAttr3pt[mTag]['gamma']

                    for z3 in dispListAvg:
                        for gamma in gammaList:
                            dkey = (mTag,z3,gamma)
                            self.bins[fit][dkey] = {}
                            self.mean[fit][dkey] = {}

                            # Need separate for-loop for these because all values of ri
                            # are needed in each ri-iteration further down
                            tOpt   = {}
                            dkeyP  = {}
                            optFit = {}
                            for ri in self.RI:
                                tOpt[ri]   = getSelectedFits(fit,ri,mTag,z3)
                                dkeyP[ri]  = getOptimalPlatKeys(tOpt[ri],z3,gamma)
                                optFit[ri] = getOptimalPlatFits(fit,ri,mTag,dkeyP[ri])

                            # The off-center values are only needed for the real part
                            for pz in self.pzPosOff:
                                mT = mTag_0 if 'p0' in pz else mTag
                                platBins['Re'][pz] = self.plat.Mbins[fit]['Re'][mT][dkeyP['Re'][pz]][optFit['Re'][pz]]

                            # Evaluate the ITDs
                            for ri in self.RI:                
                                # The 'center' value is needed for both real and imaginary
                                platBins[ri]['c'] = self.plat.Mbins[fit][ri][mTag][dkeyP[ri]['c']][optFit[ri]['c']]

                                # Still use the Real part if z3 = 0 and/or mom = 0 
                                self.bins[fit][dkey][ri] = ( (platBins[ri]  ['c']    / platBins['Re']['z0']) *
                                                             (platBins['Re']['p0z0'] / platBins['Re']['p0']) )

                                self.mean[fit][dkey][ri] = jackknife.mean(self.bins[fit][dkey][ri],
                                                                          Nbins = self.Nbins, Nspl=1)

                    print('%s ITD for momentum %s completed'%(fit,mom))
                # End for momentum
            # End for fitLabels
        # End evaluatePlateauITD() -------------


        def evaluateSummationITD(fitLabels):
            # This function is specific to the summation fits
            def getOptimalSummKeys(z3_c,gamma):
                z3_0 = 0
                kp = {}
                for pz in self.pzPos:
                    zT = z3_0 if 'z0' in pz else z3_c
                    kp[pz] = (zT,gamma)
                return kp
            #-------------------

            # Zero momentum
            mTag_0 = tags.momString([0,0,0])

            # Temporary variable
            summBins = {}
            for ri in self.RI:
                summBins[ri] = {}

            for fit in fitLabels: # These are just the summation fit labels!
                for mom in self.momAvg:
                    mTag = tags.momString(mom)
                    dispListAvg = self.dispAvg[mTag]
                    gammaList   = self.dSetAttr3pt[mTag]['gamma']

                    for z3 in dispListAvg:
                        for gamma in gammaList:
                            dkey = (mTag,z3,gamma)
                            self.bins[fit][dkey] = {}
                            self.mean[fit][dkey] = {}

                            # Need separate for-loop for these because all values of ri
                            # are needed in each ri-iteration further down
                            tOpt   = {}
                            dkeyS  = {}
                            for ri in self.RI:
                                tOpt[ri]   = getSelectedFits(fit,ri,mTag,z3)
                                dkeyS[ri]  = getOptimalSummKeys(z3,gamma)

                            # The off-center values are only needed for the real part
                            for pz in self.pzPosOff:
                                mT = mTag_0 if 'p0' in pz else mTag
                                fpT = 'M_tL%d'%(tOpt['Re'][pz]) # Tag of the matrix element at the selected fit time
                                summBins['Re'][pz] = self.summ.bins[fit][fpT]['Re'][mT][dkeyS['Re'][pz]]                               

                            # Evaluate the ITDs
                            for ri in self.RI:                
                                # The 'center' value is needed for both real and imaginary
                                fpTc = 'M_tL%d'%(tOpt[ri]['c'])
                                summBins[ri]['c'] = self.summ.bins[fit][fpTc][ri][mTag][dkeyS[ri]['c']]

                                # Still use the Real part if z3 = 0 and/or mom = 0 
                                self.bins[fit][dkey][ri] = ( (summBins[ri]  ['c']    / summBins['Re']['z0']) *
                                                             (summBins['Re']['p0z0'] / summBins['Re']['p0']) )

                                self.mean[fit][dkey][ri] = jackknife.mean(self.bins[fit][dkey][ri],
                                                                          Nbins = self.Nbins, Nspl=1)

                    print('%s ITD for momentum %s completed'%(fit,mom))
                # End for momentum
            # End for fitLabels
        # End evaluateSummationITD() -------------


        for fType in self.fitTypes.keys():
            if fType == 'Plateau':
                evaluatePlateauITD(self.fitTypes['Plateau'])
            elif fType == 'Summation':
                evaluateSummationITD(self.fitTypes['Summation'])

        print('ITD evaluation completed')
    # End evaluate() -------------
