import scipy as SP
import pylab as PL
import os
import cPickle
import scipy.misc
import sys
import pdb
import pyplot
import copy
import scipy.stats as st

from models import *
from data import *

if __name__ == '__main__':
    PL.ion()
    samples_res = 178
    samples_sus = 280
    Nsamp = 10
    #32.49 coverage over both
    Nread = 20
    chrom_len = 17645844
    #recombination rate: 2.0/chrom_len
    Rc=2.0*1.0/ (chrom_len)

    xc = 20129204
    DX = SP.linspace(0,2E6,2000)

    #probability of uneven number of recombination events
    Podd = 0.5 * (1-SP.exp(-2.0*DX*Rc))
    Peven = 1-Podd

    
    SAMPr = SP.zeros([Podd.shape[0],Nsamp])
    SAMPk = SP.zeros([Podd.shape[0],Nsamp])
    #sample a pool from these
    for i in xrange(Podd.shape[0]):
        #theoretical rate
        r = Peven[i]
        #sample pool
        binR=st.binom(samples_res,r)
        for s in xrange(Nsamp):
            rs=binR.rvs(1)
            binSeq=st.binom(Nread,SP.double(rs)/samples_res)
            kSeq = binSeq.rvs(1)
            SAMPr[i,s] = rs
            SAMPk[i,s] = kSeq

    #1. plot theoretical curve
    PL.figure()
    PL.subplot(411)
    PL.plot(DX,Peven)
    if 1:
        PL.subplot(412)
        #2. plot samples within pool
        PL.plot(DX,SAMPr/samples_res)
    if 1:
        PL.subplot(413)
        PL.plot(DX,SAMPk/Nread)

        bin_args = {'window':500E3,'step':10E3,'average':True}        
        pos = SP.concatenate((DX[:,SP.newaxis],DX[:,SP.newaxis]),axis=1)
        fres_bin,p = bin_data(SAMPk/Nread,pos,**bin_args)
        PL.subplot(414)
        PL.plot(p,fres_bin)
    
