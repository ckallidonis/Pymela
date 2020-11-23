'''
Created on Nov.20, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

This file contains definitions of conventions about the keys of the objects(dictionaries) of the JSON-format input files.
If these conventions are modified, they are also modified just once within the Pymela package
'''


# Main Object Conventions
analysisInfoTag = 'Analysis Info'
c2ptDataInfoTag = '2pt Data Info'
ensembleInfoTag = 'Ensemble Info'

# What is expected in the JSON input file, based on the type of run/test
inputInfoTags = {'2pt analysis': [analysisInfoTag, c2ptDataInfoTag, ensembleInfoTag]}

# What is expected in each object of the JSON input file
expectedInput = {c2ptDataInfoTag: ['Data Main Directory', 'Nt', 'Source Operator File', 'Sink Operator File',
                                   'Source Nrows', 'Sink Nrows', 'Compute X-rows'],
                 analysisInfoTag: ['Phasing Tag', 't0', 'Momentum', 'Nvec', 'Ncfg'],
                 ensembleInfoTag: ['Tag', 'L', 'T', 'alat fm', 'mpi GeV', 'mN GeV']}
