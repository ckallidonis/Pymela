'''
Created on Nov.24, 2020
@author: Christos Kallidonis
Copyright (C) 2020. All rights reserved.

Contains functions that define the directory and file format of input/output files
'''

def getTwoPointFileNameASCII(phM,t0Tag,src,snk,srow,mFTag,nvec):
        filePre = 'corr_2pt.baryon.n%d'%(nvec)
    
        srcTag = 'src_%s_%d'%(src,srow)
        snkTag = 'snk_%s_%d'%(snk,srow)
    
        FileName = '%s.%s.%s.%s.%s.%s.dat'%(filePre,phM,t0Tag,srcTag,snkTag,mFTag)
    
        return FileName
#--------------------------

def getTwoPointDirASCII(mainDir,t0Tag,mFTag):
        return '%s/%s/%s' % (mainDir,t0Tag,mFTag)
