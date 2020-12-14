'''
Created on Nov.20, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

This file contains functions related to I/O
'''

import json
import pymela.io.io_conventions as ioConv

def parse(ini_file):
    with open(ini_file) as f:
        data = json.load(f)
    return data
#-------------------------------

def dumpDictObject(dictObj,printMessage=''):
    if printMessage != '':
        print(printMessage)
    print(json.dumps(dictObj, indent = 2, sort_keys=True))
#-------------------------------

def makeInputChecks(runType, ioDict):
    for infoTag in ioConv.inputInfoTags[runType]:
        if infoTag not in ioDict.keys():
            raise ValueError('Expected object %s in JSON input file' % (infoTag))

        infoDict = ioDict[infoTag]

        for key in ioConv.expectedKeys[infoTag]:
            if key not in infoDict:
                raise ValueError('Expected entry "%s" in object "%s" of JSON input file' % (key,infoTag))

            if infoTag in ioConv.expectedSubKeys.keys() and key in ioConv.expectedSubKeys[infoTag].keys():
                for val in ioConv.expectedSubKeys[infoTag][key]:
                    for subDict in infoDict[key]:
                        if val not in subDict.keys():
                            raise ValueError('Expected entry "%s" in sub-object "%s/%s" of JSON input file' % (val,infoTag,key))

        if( (infoTag == ioConv.c2ptDataInfoTag or infoTag == ioConv.c3ptDataInfoTag or
             infoTag == ioConv.ratioInfoTag) and infoDict["Write HDF5 Output"]):
            if "HDF5 Output File" not in infoDict:
                raise ValueError('Got "Write HDF5 Output"=True for %s, but no file is provided. Please define "HDF5 Output File".' %(infoTag))

        
        if infoTag in ioConv.optionalKeys.keys():
            for key in ioConv.optionalKeys[infoTag]:
                if key in infoDict.keys():
                    for subDict in infoDict[key]:
                        for val in ioConv.expectedSubKeys[infoTag][key]:
                            if val not in subDict.keys():
                                raise ValueError('Expected entry "%s" in sub-object "%s/%s" of JSON input file' % (val,infoTag,key))

                            if val in ioConv.expectedSubSubKeys.keys():
                                for subVal in ioConv.expectedSubSubKeys[val]:
                                    if subVal not in subDict[val].keys():
                                        raise ValueError('Expected entry "%s" in sub-object "%s/%s/%s" of JSON input file' % (subVal,infoTag,key,val))

                        if subDict['Write HDF5 Output'] and 'HDF5 Output File' not in subDict:
                            raise ValueError('Got "Write HDF5 Output"=True for %s/%s with label %s, but no file is provided. Please define "HDF5 Output File".' %(infoTag,key,subDict['Label']))




#-------------------------------
