'''
Created on Nov.20, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

This file contains functions related to I/O
'''

import json
import io_conventions as ioConv

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

def makeInputChecks(runType,infoTag,infoDict):
    if infoTag not in ioConv.inputInfoTags[runType]:
        raise ValueError('Expected object %s in JSON input file' % (infoTag))
        
    for key in ioConv.expectedKeys[infoTag]:
        if key not in infoDict:
            raise ValueError('Expected entry "%s" in object "%s" of JSON input file' % (key,infoTag))

        if infoTag in ioConv.expectedSubKeys.keys() and key in ioConv.expectedSubKeys[infoTag].keys():
            for val in ioConv.expectedSubKeys[infoTag][key]:
                for subDict in infoDict[key]:
                    if val not in subDict.keys():
                        raise ValueError('Expected entry "%s" in sub-object "%s/%s" of JSON input file' % (val,infoTag,key))
#-------------------------------
