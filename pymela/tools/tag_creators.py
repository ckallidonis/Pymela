'''
Created on Nov.23, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Helper functions that create tags used as keys in various objects
'''

valSgnPN = lambda i: ("+" if i > 0 else "") + str(i)
valSgnN  = lambda i: ("" if i > 0 else "") + str(i)

def momH5(mom,jc = '_'):
    return 'mom'+ jc + jc.join([valSgnPN(i) for i in mom])

def momVec(mStr,sep=','):
    return [int(i) for i in (mStr.split(sep))]

def momFile(mom,jc= '.'):
    return 'momXYZ'+ jc + jc.join([str(i) for i in mom])

def momString(mom,sep=','):
    return sep.join([str(i) for i in mom])

def t0(t0):
    return 't0_%d'%(t0)

def tsep(tsep):
    return 'tsnk_%d'%(tsep)

def src_snk(opPair):
    return 'src_snk_%s_%s'%(opPair[0],opPair[1])

def row(row):
    return 'row_%d'%(row)

def disp(z3):
    dL = lambda i: ("z" if i != 0 else "") + ("+" if i > 0 else "") + str(i)
    return 'disp_%s'%(dL(z3))
