'''
Created on Nov.23, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Helper functions that create tags used as keys in various objects
'''

valSgnPN = lambda i: ("+" if i > 0 else "") + str(i)
valSgnN  = lambda i: ("" if i > 0 else "") + str(i)

def momentum(mom):
    return 'p.' + '.'.join([str(i) for i in mom])

def momFile(mom):
    return 'momXYZ.' + '.'.join([str(i) for i in mom])

def momString(mom):
    return ','.join([valSgnPN(i) for i in mom])

def t0(t0):
    return 't0_%d'%(t0)
