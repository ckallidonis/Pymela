'''
Created on Nov.20, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

This file contains definitions of conventions about the keys of the objects(dictionaries) of the JSON-format input files.
If these conventions are modified, they are also modified just once within the Pymela package
'''


# Main Object Conventions
analysisInfoTag  = 'Analysis Info'
c2ptDataInfoTag  = '2pt Data Info'
c3ptDataInfoTag  = '3pt Data Info'
ensembleInfoTag  = 'Ensemble Info'
effEnergyInfoTag = 'Effective Energy Info'
ratioInfoTag     = 'Ratio Info'
ratioFitInfoTag  = 'Ratio Fitting'
ITDInfoTag       = 'ITD'

# What is expected in the JSON input file, based on the type of run/test
inputInfoTags = {'2pt analysis': [analysisInfoTag, c2ptDataInfoTag, ensembleInfoTag],
                 '3pt analysis': [analysisInfoTag, c3ptDataInfoTag, ensembleInfoTag], 
                 'Effective Energy': [analysisInfoTag, c2ptDataInfoTag, ensembleInfoTag, effEnergyInfoTag],
                 'Compute ratio': [analysisInfoTag, c2ptDataInfoTag, c3ptDataInfoTag, ensembleInfoTag, ratioInfoTag], 
                 'Fit Ratio': [analysisInfoTag, c2ptDataInfoTag, c3ptDataInfoTag, ensembleInfoTag, ratioInfoTag, ratioFitInfoTag],
                 'Compute rITD': [analysisInfoTag, c2ptDataInfoTag, c3ptDataInfoTag, ensembleInfoTag, ratioInfoTag, ratioFitInfoTag,ITDInfoTag]
                }

# What is expected in each object of the JSON input file
expectedKeys = {c2ptDataInfoTag:  ['Data Main Directory', 'Datasets', 'Data Source', 'Write HDF5 Output'],
                c3ptDataInfoTag:  ['Data Main Directory', 'Datasets', 'Data Source', 'Write HDF5 Output'],
                analysisInfoTag:  ['Phasing Tag', 'Nvec', 'Binsize'],
                ensembleInfoTag:  ['Tag', 'L', 'T', 'alat fm', 'mpi GeV', 'mN GeV'],
                effEnergyInfoTag: ['HDF5 Output File', 'Fitting'],
                ratioInfoTag    : ['Write HDF5 Output'],
                ratioFitInfoTag : [],
                ITDInfoTag: ['Optimal Fits']}

expectedSubKeys = {c2ptDataInfoTag: {'Datasets': ['Mom List', 'Ncfg', 't0','Nt','Interpolating Operators File',
                                                  'Nrows','Compute X-rows']
                                    },
                   c3ptDataInfoTag: {'Datasets': ['Mom List', 'Ncfg', 't0','tsep','disp','Interpolating Operators File',
                                                  'Insertion Operators', 'Nrows','Compute X-rows']
                                    },                 
                   effEnergyInfoTag: {'Fitting': ['Type', 'Ranges']},
                   ratioFitInfoTag: {'Plateau'  : ['Type','Label','Chi Criterion','Write HDF5 Output'],
                                     'Summation': ['Type','Label','tsepLow', 'Fit Bands']}
                  }

# Not optimal, perhaps, but a quick and dirty solution
expectedSubSubKeys = {'Fit Bands': ['Evaluate', 'Npoints']}

optionalKeys = {ratioFitInfoTag : ['Plateau','Summation','Two-state']}