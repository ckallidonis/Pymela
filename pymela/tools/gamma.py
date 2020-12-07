'''
Created on Dec. 7, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

This file contains information related to the gamma insertion currents
'''

def insertionMap(gamma):

    gMatMap =  {'gt'  : ('b_b0xDA__J0_A1pP', 1),
               'gxg5': ('a_a1xDA__J1_T1pM', 1),
               'gyg5': ('a_a1xDA__J1_T1pM', 3),
               'gzg5': ('a_a1xDA__J1_T1pM', 2),
               'gtg5': ('pion_pion_2xDA__J0_A1mM', 1),
               'gxgy': ('b_b1xDA__J1_T1pP', 2),
               'gxgz': ('b_b1xDA__J1_T1pP', 3),
               'gxgt': ('rho_rho_2xDA__J1_T1mP', 3),
               'gygz': ('b_b1xDA__J1_T1pP', 1),
               'gygt': ('rho_rho_2xDA__J1_T1mP', 1),
               'gzgt': ('rho_rho_2xDA__J1_T1mP', 2)}

    return gMatMap[gamma]
#--------------------------------------------------------