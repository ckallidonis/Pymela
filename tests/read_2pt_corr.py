'''
Created on Nov.20, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Simple test that read two-point correlation functions
'''

import sys, os
import optparse

# Add package path to sys.path. This allows us to run this tests script from any directory, without import issues
file_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(file_path+'/../')
fileName = __file__.split('/')[-1]

# Import local modules
import pymela.io.json_io as JSONio
import pymela.io.io_conventions as ioConv
from pymela.twopointcorr import TwoPointCorrelator

runType = '2pt analysis'

# Avoid writing the compiled files
sys.dont_write_bytecode = True


# Parse command line options
usage = "usage: %prog [options] "
opt_parser = optparse.OptionParser(usage)

opt_parser.add_option("-i", "--input_file", type="string", default='',
                  help='Input file in JSON format')

(options, args) = opt_parser.parse_args()

input_file=options.input_file
if input_file == '':
    raise ValueError('--input_file option must be set')


# Parse and print out the input file
ioDict = JSONio.parse(input_file)
JSONio.dumpDictObject(ioDict,'\n%s - Got the following Input:' %(fileName))

# Make cheks on Input data
JSONio.makeInputChecks(runType, ioDict)


c2pt_dataInfo = ioDict[ioConv.c2ptDataInfoTag]
analysisInfo = ioDict[ioConv.analysisInfoTag]
ensembleInfo = ioDict[ioConv.ensembleInfoTag]



c2pt = TwoPointCorrelator(dataInfo = c2pt_dataInfo, analysisInfo = analysisInfo)
c2pt.printInfo()

c2pt.getData()

# Perform Statistical / Jackknife Analysis
c2pt.doStatistics()

# Write the output in HDF5 format
if c2pt_dataInfo['Write HDF5 Output']:
    c2pt.writeHDF5()
