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
    #1. parse command line options
    if 1:
      #qual 1,000
      sys.argv = ['', '-f','./../data/vcf-snp-files/chr07_res_sus_pass_qual_1000.recode.vcf','-o','./out','-c','7','--star', 0.0E7, '--stop', 10E7, 'recalc']

      #preprocessing parameters
      preprocess_params = {'max_coverage':100, 'filter_dash':True, 'filterN':True, 'enforce_match_major': True, 'enforce_match_minor': False,'aligned':True}

    options   = parse_options(sys.argv)
    data_file = options.data_file
    out_file  = options.out_file
    #data file and caching file:
    data_file_cached = data_file +'.pickle'

    samples_res = 178
    samples_sus = 280

    #causal locus
    xc = 20129204

    preprocess_cache_file = './out/preprocess.pickle'
    #1. load raw data
    fileName, fileExtension = os.path.splitext(data_file)
    if fileExtension=='.vcf':
        data = load_data_vcf(data_file,data_file_cached)
    else:
        data = load_data(data_file,data_file_cached)

    recalc = True
    if not os.path.exists(preprocess_cache_file) or recalc:
        preprocess_data(data,**preprocess_params)
        cPickle.dump({'data':data},open(preprocess_cache_file,'wb'))
    else:
        R = cPickle.load(open(preprocess_cache_file,'rb'))
        globals().update(R)

    #4. filter data
    Iok = SP.ones(data['pos'].shape[0],dtype='bool')
    if options.chrom!='-':
      Iok = Iok & (data['pos'][:,0]==int(options.chrom))
    if options.start!='-':
      Iok = Iok & (data['pos'][:,1]>=int(options.start))
    if options.stop!='-':
      Iok = Iok & (data['pos'][:,1]<int(options.stop))

    #filter data by genomic position:
    filter_data(data,Iok)
    pos = data['pos']

    if 0:
        pdb.set_trace()
        Lres = eval_binomial(data,index_res=0,p=[0.5],fields=['counts_res'])
        Lsus = eval_binomial(data,index_res=1,p=[0.5],fields=['counts_sus'])

        Lres_bin,pos_bin1 = bin_data(Lres,pos,w=100E3)
        Lsus_bin,pos_bin2 = bin_data(Lsus,pos,w=100E3)

        PL.plot(pos_bin1,Lres_bin,'r.-')
        PL.plot(pos_bin2,Lsus_bin,'r.-') 
        

    if 1:
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

        if 0:
            Iok = L0>(-40)
            Lres = Lres[Iok]
            Lsus = Lsus[Iok]
            L0   = L0[Iok]
            pos  = pos[Iok]
        
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


    if 0:
        Nhom_res,p = bin_data(1.0*(data['counts_res'][:,1]==0),pos,w=100E3)
        Nhom_sus,p = bin_data(1.0*(data['counts_sus'][:,1]==0),pos,w=100E3)
        Nhet_res,p = bin_data(1.0*(data['counts_res'][:,1]!=0),pos,w=100E3)
        Nhet_sus,p = bin_data(1.0*(data['counts_sus'][:,1]!=0),pos,w=100E3)

        R_res = Nhet_res/(Nhom_res+Nhet_res)
        R_sus = Nhet_sus/(Nhom_sus+Nhet_sus)
        
        p0=PL.plot(p,R_res,'r.-')
        p1=PL.plot(p,R_sus,'b.-')
        p2=PL.plot(p,R_res-R_sus,'k.-')
        PL.legend([p0,p1,p2],['res','sus','delta'])
        

    if 0:
        k_res = data['counts_res'][:,0]
        N_res = data['counts_res'].sum(axis=1)
        k_sus = data['counts_sus'][:,0]
        N_sus = data['counts_sus'].sum(axis=1)

        pos_bin = pos[:,1]
        Irel = (SP.absolute(pos_bin-xc)<0.5E6)
        m_res = k_res[Irel]/N_res[Irel]
        m_sus = k_sus[Irel]/N_sus[Irel]

        PL.figure()
        p0=PL.plot(pos_bin[Irel],m_res,'b.-')
        p1=PL.plot(pos_bin[Irel],m_sus,'r.-')
        p2= PL.plot(pos_bin[Irel],m_sus-m_res,'k.-')
        PL.legend([p0,p1,p2],['res','sus','delta'])


    if 0:       
        k_res = data['counts_res'][:,0]
        N_res = data['counts_res'].sum(axis=1)
        k_sus = data['counts_sus'][:,0]
        N_sus = data['counts_sus'].sum(axis=1)

        k_res,pos_bin = bin_data(k_res,pos,w=70E3)
        N_res,pos_bin = bin_data(N_res,pos,w=70E3)
        N_sus,pos_bin = bin_data(N_sus,pos,w=70E3)
        k_sus,pos_bin = bin_data(k_sus,pos,w=70E3)
        
        Irel = (SP.absolute(pos_bin-xc)<1E10)
        m_res = k_res[Irel]/N_res[Irel]
        m_sus = k_sus[Irel]/N_sus[Irel]

        PL.figure()
        PL.plot(pos_bin[Irel],m_res,'b.')
        PL.plot(pos_bin[Irel],m_sus,'r.')
      

    
