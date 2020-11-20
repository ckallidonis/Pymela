'''
Created on Nov.20, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Class definition that holds data and manipulates Two-point correlation functions
'''

import io.json_io as JSONio

class TwoPointCorrelator():
    def __init__(self, dataInfo, analysisInfo):
        self.dataInfo = dataInfo
        self.analysisInfo = analysisInfo

    def printInfo(self):
        JSONio.dumpDictObject(self.dataInfo, '\nTwoPointCorrelator - Got the following 2pt data Info:')
