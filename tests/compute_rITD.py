'''
Created on Nov.20, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Test that evaluates reduced Ioffe-time distributions from two- and three-point functions
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
from pymela.threepointcorr import ThreePointCorrelator
from pymela.ratio import ThreeToTwoPointCorrRatio
from pymela.plateau_fit import PlateauFit
from pymela.summation_fit import SummationFit

runType = 'Compute rITD'

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
c3pt_dataInfo = ioDict[ioConv.c3ptDataInfoTag]
analysisInfo = ioDict[ioConv.analysisInfoTag]
ensembleInfo = ioDict[ioConv.ensembleInfoTag]
ratioInfo = ioDict[ioConv.ratioInfoTag]
ratioFitInfo = ioDict[ioConv.ratioFitInfoTag]

# Read the two-point functions, perform statistical/Jackknife analysis
c2pt = TwoPointCorrelator(dataInfo = c2pt_dataInfo, analysisInfo = analysisInfo)
c2pt.printInfo()
c2pt.getData()
c2pt.doStatistics()

# Write the output in HDF5 format
if c2pt_dataInfo['Write HDF5 Output']:
    c2pt.writeHDF5()
#------------------------------------------------

# Read the three-point functions, perform statistical/Jackknife analysis
c3pt = ThreePointCorrelator(dataInfo = c3pt_dataInfo, analysisInfo = analysisInfo)
c3pt.printInfo()
c3pt.getData()
c3pt.doStatistics()

# Write the output in HDF5 format
if c3pt_dataInfo['Write HDF5 Output']:
   c3pt.writeHDF5()
#------------------------------------------------

# Define and evaluate the three- to two-point function ratios
ratio = ThreeToTwoPointCorrRatio(c2pt = c2pt, c3pt = c3pt, dataInfo = ratioInfo, analysisInfo = analysisInfo)
ratio.evaluate()

# Write the output in HDF5 format
if ratioInfo['Write HDF5 Output']:
   ratio.writeHDF5()
#------------------------------------------------

# Perform Plateau fits on the plain ratio
# if 'Plateau' in ratioFitInfo:
#     print('Will perform Plateau Fits on the Plain ratio')
#     plat = PlateauFit(ratio=ratio, ratioType='plain', fitInfo = ratioFitInfo['Plateau'], analysisInfo = analysisInfo)
#     plat.performFits()
#     plat.writeHDF5()

# Perform Linear fits on the summed ratio
if 'Summation' in ratioFitInfo:
    print('Will perform Fits on the summed ratio')
    summ = SummationFit(ratio=ratio, ratioType='sum', fitInfo = ratioFitInfo['Summation'], analysisInfo = analysisInfo)
#    summ.performFits()
#    plat.writeHDF5()

