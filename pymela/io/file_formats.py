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

def getTwoPointDirASCII(mainDir,phDir,t0Tag,mFTag):
        return '%s/%s/%s/%s' % (mainDir,phDir,t0Tag,mFTag)
#--------------------------

def getThreePointDirASCII(mainDir,phDir,t0Tag,tsepTag,mFTag):
    return '%s/%s/%s/%s/%s' % (mainDir,phDir,t0Tag,mFTag,tsepTag)
#--------------------------

def getThreePointFileNameASCII(phM,t0Tag,tsepTag,srcOp,snkOp,row,insOp,insRow,mFTag,dispTag,nvec):
    filePre = 'corr_3pt.baryon.n%d'%(nvec)

    srcTag = 'src_%s_%d'%(srcOp,row)
    snkTag = 'snk_%s_%d'%(snkOp,row)
    insTag = 'ins_%s_%d'%(insOp,insRow)

    FileName = '%s.%s.%s.%s.%s.%s.%s.%s.%s.dat'%(filePre,phM,t0Tag,tsepTag,srcTag,snkTag,insTag,dispTag,mFTag)

    return FileName
#--------------------------
