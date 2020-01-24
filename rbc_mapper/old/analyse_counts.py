import scipy as SP
import pylab as PL
import os
import cPickle
import scipy.misc
import sys
import pdb
from scipy.special import gammaln, betaln

def lnchoose(n, m):
  nf = gammaln(n + 1)
  mf = gammaln(m + 1)
  nmmnf = gammaln(n - m + 1)
  return nf - (mf + nmmnf)

def log_binomial(k,N,p,norm=False):
    if norm:
        lZ = lnchoose(N,k)
    else:
        lZ = 0
    pdf = lZ + k * SP.log(p) + (N-k)*SP.log(1-p)
    return pdf

def trivial_phasing(counts_res,counts_sus,do_phase=True):
    #separate in N and k
    #k is the major allele in res!
    iN = SP.arange(counts_res.shape[0])
    if do_phase:
        ik = counts_res.argmax(axis=1)
    else:
        ik = SP.ones(counts_res.shape[0],dtype='int')
        

    #get total count and k for res and sus
    N_res = counts_res.sum(axis=1)
    k_res = counts_res[iN,ik]
    N_sus = counts_sus.sum(axis=1)
    k_sus = counts_sus[iN,ik]
    return N_res,k_res,N_sus,k_sus



def binomial_test_stepwise(N_res,k_res,N_sus,k_sus,pos,p_res=0.99,p_sus=0.33,p0=0,p1=SP.inf,step_size=10E4,d_max=0,norm=False):
    #1. generate total null distribution
    lpdf_res0 = log_binomial(k_res,N_res,0.5,norm=norm)
    lpdf_sus0 = log_binomial(k_sus,N_sus,0.5,norm=norm)
    L0 = (lpdf_res0+lpdf_sus0).sum()

    p0 = max(p0,pos.min())
    p1 = min(p1,pos.max())
    #d_max = (p1-p0)
    #2. test loci
    x = p0
    R = []
    P = []
    Z = []
    while x<p1:
        x+= step_size
        d = SP.absolute(pos-x)
        P_res = 1-d/d_max*0.5
        P_res[P_res<0.5] = 0.5
        P_sus = 0.33 + d/d_max*(0.5-0.33)
        P_sus[P_sus>0.5] = 0.5
        S = log_binomial(k_res,N_res,P_res,norm=norm) + log_binomial(k_sus,N_sus,P_sus,norm=norm)
        #calc posterior
        Sp = S.sum()
        nom   = SP.exp(0)
        denom = SP.exp(0) + SP.exp(L0-Sp)
        post = nom/denom      
        if 0:
            pdb.set_trace()
            PL.figure(1),PL.clf(),PL.plot(P_sus),PL.plot(P_res)
            PL.figure(2),PL.clf(),PL.plot(S)
            pdb.set_trace()
        R.append(Sp-L0)
        P.append(x)
        Z.append(post)
    R = SP.array(R)
    P = SP.array(P)

    return P,R,Z
        
    

def binomial_test(N_res,k_res,N_sus,k_sus,p_res=0.99,p_sus=0.33):
    #binomail test on position-wide scale
    lpdf_res = log_binomial(k_res,N_res,p_res)
    lpdf_sus = log_binomial(k_sus,N_sus,p_sus)

    #bacground models
    lpdf_res0 = log_binomial(k_res,N_res,0.5)
    lpdf_sus0 = log_binomial(k_sus,N_sus,0.5)

    score_res = lpdf_res-lpdf_res0
    score_sus = lpdf_sus - lpdf_sus0
    score = score_res + score_sus

    return score,score_res,score_sus



def bin_data(d,pos,w=400E3,step=10E3):
    # 1. binning
    d_bin = []
    pos_bin = []
    p0 = 0
    pmax = pos[:,1].max()
    while p0 < pmax:
        p1 = p0 + w
        Ib = (p0<=pos[:,1]) & (pos[:,1]<p1)
        d_bin.append(d[Ib,:].sum(axis=0))
        pos_bin.append(pos[Ib,1][0])
        p0 += step
    d_bin = SP.array(d_bin)
    pos_bin    = SP.array(pos_bin)
    return d_bin,pos_bin
             

if __name__ == '__main__':
    data_file = './../data/all_alleles.txt'
    data_file_cached = data_file +'.pickle'

    max_coverage = 300

    recalc = 'recalc' in sys.argv
    if not os.path.exists(data_file_cached) or recalc:
        M = SP.loadtxt(data_file,dtype='str')

        i = 0
        pos = SP.asarray(M[1::,i:i+2],dtype='int')
        i+=2
        alleles_sus = M[1::,i:i+2]
        i+=2
        counts_sus  = SP.asarray(M[1::,i:i+2],dtype='float')
        i+=2
        alleles_res = M[1::,i:i+2]
        i+=2
        counts_res  =SP.asarray(M[1::,i:i+2],dtype='float')
        i+=2
        #skip both
        i+=4
        #parents
        #resistant pool
        alleles_UM = M[1::,i:i+2]
        i+=2
        counts_UM = SP.asarray(M[1::,i:i+2],dtype='float')
        i+=2
        alleles_TN1 = M[1::,i:i+2]
        i+=2
        counts_TN1 = SP.asarray(M[1::,i:i+2],dtype='float')

        #combined allesl  in pool
        alleles_pool = SP.concatenate((alleles_sus,alleles_res),axis=1)
        counts_pool  = SP.concatenate((counts_sus,counts_res),axis=1)
        counts_parent= SP.concatenate((counts_UM,counts_TN1),axis=1)
        
        Iok = SP.ones((alleles_pool.shape[0]),dtype='bool')       
        Iok = Iok & (counts_pool.sum(axis=1)<max_coverage)
        #check homocygocity of parents
        Iok = Iok & (alleles_UM[:,1]=='N') & (alleles_TN1[:,1]=='N')
        Iok = Iok & (counts_parent.sum(axis=1)>8)
        Iok = Iok & (counts_pool.sum(axis=1)>8)

        alleles_UM = alleles_UM[Iok]
        alleles_TN1= alleles_TN1[Iok]
        counts_UM = counts_UM[Iok]
        counts_TN1= counts_TN1[Iok]
        alleles_res=alleles_res[Iok]
        alleles_sus=alleles_sus[Iok]
        counts_sus = counts_sus[Iok]
        counts_res = counts_res[Iok]
        pos = pos[Iok]
        

        #now align based [non-res,res]
        aUM = alleles_UM[:,0]
        aTN1= alleles_TN1[:,0]

        Iflip_res = (alleles_res[:,1]==aTN1)
        Iflip_sus = (alleles_sus[:,1]==aTN1)

        Iflip_res = Iflip_res | (alleles_res[:,0]==aUM)
        Iflip_sus = Iflip_sus | (alleles_sus[:,0]==aUM)

        if 1:
            alleles_sus[Iflip_sus] = alleles_sus[Iflip_sus][:,[1,0]]
            counts_sus[Iflip_sus] = counts_sus[Iflip_sus][:,[1,0]]

            alleles_res[Iflip_res] = alleles_res[Iflip_res][:,[1,0]]
            counts_res[Iflip_res] = counts_res[Iflip_res][:,[1,0]]

        
        RV = {'pos':pos,'counts_res':counts_res,'counts_sus':counts_sus,'alleles_sus':alleles_sus,'alleles_res':alleles_res}
        cPickle.dump(RV,open(data_file_cached,'wb'),-1)
    else:
        RV = cPickle.load(open(data_file_cached,'rb'))
        globals().update(RV)


    if 0:
        #plotting of sus/res and things
        bin_sus,pos_bin = bin_data(counts_sus,pos,w=400E3)
        bin_res,pos_bin = bin_data(counts_res,pos,w=400E3)

        N_res,k_res,N_sus,k_sus = trivial_phasing(bin_res,bin_sus,do_phase=True)

        PL.figure()
        PL.plot(pos_bin,bin_res[:,1]/bin_res.sum(axis=1))
        PL.plot(pos_bin,bin_sus[:,1]/bin_sus.sum(axis=1))
        PL.figure()
        PL.plot(pos_bin,k_res/N_res)
        PL.plot(pos_bin,k_sus/N_sus)

    

    if 1:
        #global model with changing rates from recombination
        p0 = 0.0E7
        p1 = 4E7

#        p0 = 1.9E7
#        p1 = 2.5E7

        Irel = (pos[:,1]>=p0) & (pos[:,1]<p1)
        pos = pos[Irel]
        counts_res = counts_res[Irel]
        counts_sus = counts_sus[Irel]
        N_res,k_res,N_sus,k_sus = trivial_phasing(counts_res,counts_sus,do_phase=True)

        rate = 1./5
        P,R,Z = binomial_test_stepwise(N_res,k_res,N_sus,k_sus,pos[:,1],d_max= rate*(pos[:,1].max()-pos[:,1].min()),step_size=1E4)

        PL.plot(P,R,'k.')
        PL.savefig('recombination_model.png')

    if 0:
        bin_sus,pos_bin = bin_data(counts_sus,pos,w=400E3)
        bin_res,pos_bin = bin_data(counts_res,pos,w=400E3)

        PL.figure(figsize=[10,6])
        PL.subplot(211)
        PL.plot(pos_bin,bin_sus[:,1]/bin_sus.sum(axis=1))
        PL.subplot(212)
        PL.plot(pos_bin,bin_res[:,1]/bin_res.sum(axis=1))

    if 0:
        #get trivially phased data
        N_res,k_res,N_sus,k_sus = trivial_phasing(counts_res,counts_sus,do_phase=True)

        #carry out independeng binomial test
        score,score_res,score_sus = binomial_test(N_res,k_res,N_sus,k_sus)

        #plot
        PL.figure(figsize=[8,10])
        PL.title('raw scores per SNP')
        PL.subplot(311)
        PL.title('score res')
        PL.plot(pos[:,1],score_res,'k.')
        PL.subplot(312)
        PL.title('score sus')
        PL.plot(pos[:,1],score_sus,'k.')
        PL.subplot(313)
        PL.title('score both')
        PL.plot(pos[:,1],score,'k.')
        PL.savefig('raw_scores.png')

        #carry out simple binnin
        bin_scores,pos_bin = bin_data(score,pos)
        PL.figure(figsize=[8,4])
        PL.title('Binned score')
        PL.plot(pos_bin,bin_scores,'k.')
        PL.savefig('binned_scores.png')

        PL.xlim([1.9E7,2.3E7])
        PL.savefig('binned_scores_zoom.png')
    

