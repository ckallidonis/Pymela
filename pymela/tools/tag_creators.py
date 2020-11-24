'''
Created on Nov.23, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Helper functions that create tags used as keys in various objects
'''

valSgnPN = lambda i: ("+" if i > 0 else "") + str(i)
valSgnN  = lambda i: ("" if i > 0 else "") + str(i)

def momentum(zMom):
    return 'pz.%s'%(valSgnPN(zMom))

