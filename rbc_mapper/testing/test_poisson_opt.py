import scipy as SP
import pylab as PL
import os
import cPickle
import scipy.misc
import sys
import pdb
import pyplot
import copy 

from models import *
from data import *


if __name__ == '__main__':
    #causal locus
    xc = 20129204
    
    #1. parse command line options
    if 1:
      #qual 1,000
      sys.argv = ['', '-f','./../data/snp_files/oli_input_2_samples.vcf','-o','./out','-c','07','--star', xc-100E6, '--stop', xc+100E6]
      sys.argv = ['', '-f','./../data/snp_files/oli_input_4_samples.vcf','-o','./out','-c','07','--star', xc-100E6, '--stop', xc+100E6,'--with_parents']

      #preprocessing parameters
      preprocess_params = {'max_coverage':100, 'filter_dash':True, 'filterN':True, 'enforce_match_major': True, 'enforce_match_minor': False,'filter_parent':True}

    options   = parse_options(sys.argv)
    data_file = options.data_file
    out_file  = options.out_file
    #data file and caching file:
    data_file_cached = data_file +'.pickle'

    samples_res = 178
    samples_sus = 280


    recalc = False
    #1. load raw data
    fileName, fileExtension = os.path.splitext(data_file)
    if not os.path.exists(data_file_cached) or recalc:
        data = load_data_vcf(data_file,with_parents=options.with_parents)
        cPickle.dump(data,open(data_file_cached,'wb'),-1)
    else:
        data = cPickle.load(open(data_file_cached,'rb'))

    preprocess_data(data,**preprocess_params)    
    #4. filter data
    Iok = SP.ones(data['pos'].shape[0],dtype='bool')
    if options.chrom!='-':
      Iok = Iok & (data['chrom']==options.chrom)
    filter_data(data,Iok)
    chrom_len = data['pos'].max()-data['pos'].min()
    Iok = SP.ones(data['pos'].shape[0],dtype='bool')
    if options.start!='-':
      Iok = Iok & (data['pos']>=int(options.start))
    if options.stop!='-':
      Iok = Iok & (data['pos']<int(options.stop))
    filter_data(data,Iok)
    pos = data['pos']

    #get library size correction
    [LSres,LSsus] = lib_size_factors(data)
    #chromosome length:

    pos = data['pos']
    #recomb. rate for a single samples per bp
    Rc = 2.0/chrom_len
    #Distance until which 50% chance of a recombination
    D50 = 0.5/Rc

    if 1:
        if 0:
            p0 = xc-4E6
            p1 = xc+4E6
        if 1:
            p0 = pos.min()
            p1 = pos.max()
        #### new recombination base model ###
        rm = RcombinationMapping()
        rm.setCountData(data['counts_res'],data['counts_sus'],data['pos'])
        rm.setPoolSizes(samples_res,samples_sus)
        rm.setRecombination(100*1.0/chrom_len)
        [P,S,S0]=rm.score(p0,p1,step_size=50E3,window_size=500E3)
        pdb.set_trace()
        PL.figure()
        PL.subplot(211)
        PL.plot(P,S0,'k--')
        PL.plot(P,S)
        PL.subplot(212)
        PL.plot(P,S.max(axis=1)-S0)
        #PL.plot(P,S0)


    if 0:
        PL.ion()
        data['counts_both'] = LSres*data['counts_res']+LSsus*data['counts_sus']
        globals().update(data)

        bin_args = {'window':400E3,'step':10E3}
        #bin_args = {'N':10}
        k_res = counts_res[:,0]
        N_res = counts_res.sum(axis=1)
        k_sus = counts_sus[:,0]
        N_sus = counts_sus.sum(axis=1)
        k_both = counts_both[:,0]
        N_both = counts_both.sum(axis=1)
        counts_res_bin,p = bin_data(counts_res,pos,**bin_args)
        counts_sus_bin,p = bin_data(counts_sus,pos,**bin_args)
        data['counts_res_bin'] = counts_res_bin
        data['counts_sus_bin'] = counts_sus_bin

        #binomial ealuation for 1.0 and 0.33
        L0 = eval_binomial(data,index_res=0,p=[0.5],fields=['counts_both'])
        L0_bin,p = bin_data(L0,pos,**bin_args)

        eps = 1E-1
        Lres=eval_binomial(data,index_res=0,fields=['counts_res'],p=[1-eps])
        Lsus=eval_binomial(data,index_res=0,fields=['counts_sus'],p=[0.3])
        Lres0=eval_binomial(data,index_res=0,fields=['counts_res'],p=[0.5])
        Lsus0=eval_binomial(data,index_res=0,fields=['counts_sus'],p=[0.5])

        Lbal=eval_binomial(data,index_res=0,fields=['counts_both'],p=[0.5])

        Lbal_bin,p = bin_data(Lbal,pos,**bin_args)
        Lres_bin,p = bin_data(Lres,pos,**bin_args)
        Lsus_bin,p = bin_data(Lsus,pos,**bin_args)
        Lres0_bin,p = bin_data(Lres0,pos,**bin_args)
        Lsus0_bin,p = bin_data(Lsus0,pos,**bin_args)


        if 1:
            PL.subplot(211)
            pres=PL.plot(p,Lres_bin-Lres0_bin)
            psus=PL.plot(p,Lsus_bin-Lsus0_bin)
            PL.legend([pres,psus],['res','sus'])
            PL.ylim([0,PL.ylim()[1]])
            PL.subplot(212)
            pjoint =PL.plot(p,Lres_bin-Lres0_bin+Lsus_bin-Lsus0_bin)
            PL.ylim([0,PL.ylim()[1]])
            PL.legend([pjoint],['joint score'])


    if 0:
        PL.figure(1)
        Nhet_res_bin+=1
        Nhet_sus_bin+=1
        Nhom_res_bin+=1
        Nhom_sus_bin+=1
        PL.subplot(211)
        fhom_res =Nhom_res_bin/(Nhom_res_bin+Nhet_res_bin)
        fhom_sus =Nhom_sus_bin/(Nhom_sus_bin+Nhet_sus_bin) 
        #fhom_res =Nhom_res_bin/(Nhet_res_bin)
        #fhom_sus =Nhom_sus_bin/(Nhet_sus_bin) 

        PL.plot(p,fhom_res)
        PL.plot(p,fhom_sus)
        PL.subplot(212)
        PL.plot(p,fhom_res-fhom_sus)

        


    if 0:
        PL.figure(2)
        PL.plot(p,f_res_bin)
        PL.plot(p,f_sus_bin)
        PL.plot(p,f_both_bin)

        PL.figure(3,figsize=[10,10])
        PL.subplot(411)
        PL.plot(p,Lres_bin)
        PL.plot(p,Lsus_bin)
        PL.subplot(412)
        PL.plot(p,Lres0_bin)
        PL.plot(p,Lsus0_bin)
        PL.subplot(413)
        PL.plot(p,Lres_bin-Lres0_bin)
        PL.plot(p,Lsus_bin-Lsus0_bin)
        PL.subplot(414)
        PL.plot(p,Lbal_bin)
        #PL.plot(p,Lsus_bin-Lsus0_bin)

        PL.show()
        #PL.plot(p,Lsus_bin)

        


    if 0:

        PL.figure()
        PL.plot(p,f_res_bin)
        PL.plot(p,f_sus_bin)

        PL.figure()
        PL.plot(pos[:,1],f_res)
        PL.plot(pos[:,1],f_sus)

        
    if 0:
        #current model and optimization-based attempt
        data['counts_both'] = data['counts_res']+data['counts_sus']
        k_res = data['counts_res'][:,0].sum()
        N_res = data['counts_res'].sum()
        k_sus = data['counts_sus'][:,0].sum()
        N_sus = data['counts_sus'].sum()

        R_res = k_res/N_res
        R_sus = k_sus/N_sus

        #fehlermodell
        L0 = eval_binomial(data,index_res=0,p=[0.5],fields=['counts_both'])

        Lres0 = eval_binomial(data,index_res=0,p=[R_res],fields=['counts_res'])
        Lsus0 = eval_binomial(data,index_res=0,p=[R_sus],fields=['counts_sus'])

        Lres1 = eval_binomial(data,index_res=1,p=[R_res],fields=['counts_res'])
        Lsus1 = eval_binomial(data,index_res=1,p=[R_sus],fields=['counts_sus'])

        Lres = SP.concatenate((Lres0[:,SP.newaxis],Lres1[:,SP.newaxis]),axis=1).max(axis=1)
        Lsus = SP.concatenate((Lsus0[:,SP.newaxis],Lsus1[:,SP.newaxis]),axis=1).max(axis=1)

        Lres_bin,p = bin_data(Lres,pos,w=1.0*70E3)
        Lsus_bin,p = bin_data(Lsus,pos,w=1.0*70E3)
        L0_bin,p   = bin_data(L0,pos,w=1.0*70E3)
        score = Lres_bin-Lsus_bin

        if 1:
            Iok = (L0_bin>-6000)
            Lres_bin = Lres_bin[Iok]
            Lsus_bin = Lsus_bin[Iok]
            L0_bin   = L0_bin[Iok]
            score = score[Iok]
            p  = p[Iok]
        
        
        PL.figure(figsize=[12,6])
        PL.subplot(211)
        PL.title('Peak score')
        PL.plot(p,score,'b.-')
        PL.ylim([0,5000])
        PL.subplot(212)
        PL.title('Quality score (sum=0.5)')
        PL.plot(p,L0_bin,'b.-')
        PL.ylim([-25000,0])
        fn = os.path.join('./../results/figures/',os.path.basename(options.data_file)+'_maxcov%d_chrom%s' % (preprocess_params['max_coverage'],options.chrom))
        fn_figure = fn+'.png'
        fn_csv = fn+'.csv'
        PL.savefig(fn_figure)
        M = SP.concatenate((p[:,SP.newaxis],score[:,SP.newaxis],L0_bin[:,SP.newaxis]),axis=1)
        SP.savetxt(fn_csv,M)

