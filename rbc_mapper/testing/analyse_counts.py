import scipy as SP
import pylab as PL
import os
import cPickle
import scipy.misc
import sys
import pdb
import pyplot

from models import *


if __name__ == '__main__':
    #1. parse command line options
    if 1:
      sys.argv = ['', '-f','./../data/snp_files/corenonrep_alleles_sus_res_join_102711-preprocessed','-o','./out','-c','7','--star', 1.0E7, '--stop', 4E7]
      #max coverage:
      preprocess_params = {'max_coverage':100, 'filter_dash':False, 'filterN':True, 'enforce_match_major': True, 'enforce_match_minor': False,'aligned':True}

    #TODO move as command line argument
    options   = parse_options(sys.argv)
    data_file = options.data_file
    out_file  = options.out_file
    #data file and caching file:
    data_file_cached = data_file +'.pickle'

    samples_res = 178
    samples_sus = 280

    samples_res = 458
    samples_sus = 458


    #causal locus
    xc = 20129204

    #1. load raw data
    data = load_data(data_file,data_file_cached)

    #2. preprocess and filter:
    preprocess_data(data,**preprocess_params)
    
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
      
    #3. trivial phasing:
    N_res,k_res,N_sus,k_sus = trivial_phasing(data,do_phase=True)
    pos = data['pos']

    if 0:
      score,score_res,score_sus = binomial_test(N_res,k_res,N_sus,k_sus)
      Irel = (SP.absolute(pos[:,1]-xc)<2E6)
      m_res = k_res[Irel]/N_res[Irel]
      m_sus = k_sus[Irel]/N_sus[Irel]

      PL.figure()
      PL.plot(pos[Irel,1],m_res,'b.')
      PL.plot(pos[Irel,1],m_sus,'r.')
      

    if 0:
      counts_res = data['counts_res']
      counts_sus = data['counts_sus']
      #plotting of sus/res and things
      bin_sus,pos_bin = bin_data(counts_sus,pos,w=77493,step=20E3)
      bin_res,pos_bin = bin_data(counts_res,pos,w=77493,step=20E3)

      N_res,k_res,N_sus,k_sus = trivial_phasing({'counts_res':bin_res,'counts_sus':bin_sus},do_phase=False)

      PL.figure()
      PL.plot(pos_bin,bin_res[:,0]/bin_res.sum(axis=1))
      PL.plot(pos_bin,bin_sus[:,0]/bin_sus.sum(axis=1))
      PL.figure()
      PL.plot(pos_bin,k_res/N_res)
      PL.plot(pos_bin,k_sus/N_sus)

    if 1:
      counts_res = data['counts_res']
      counts_sus = data['counts_sus']

      score,score_res,score_sus = binomial_test(N_res,k_res,N_sus,k_sus)
      #calc window size of step-wise model
      chrom_len = pos[:,1].max()-pos[:,1].min()
      Nrecomb   = 2*(samples_res + samples_sus)
      w = chrom_len/Nrecomb
      bin_sus,pos_bin = bin_data(counts_sus,pos,w=w)
      bin_res,pos_bin = bin_data(counts_res,pos,w=w)

      score_sus,pos_bin = bin_data(score_sus,pos,w=w)
      score_res,pos_bin = bin_data(score_res,pos,w=w)
      score,pos_bin = bin_data(score,pos,w=w)           
      pdb.set_trace()  

    if 0:
        #4. run binomial test stepwside model:
        rate = 1./2
        P,R,Z = binomial_test_stepwise(N_res,k_res,N_sus,k_sus,pos[:,1],d_max= rate*(pos[:,1].max()-pos[:,1].min()),step_size=1E5)
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
      if 0:
        counts_res = data['counts_res']
        counts_sus = data['counts_sus']
        #plotting of sus/res and things
        bin_sus,pos_bin = bin_data(counts_sus,pos,w=10E4)
        bin_res,pos_bin = bin_data(counts_res,pos,w=10E4)
        N_res,k_res,N_sus,k_sus = trivial_phasing({'counts_res':bin_res,'counts_sus':bin_sus},do_phase=True)
      else:
        pos_bin = pos[:,1]       
      #carry out independeng binomial test
      score,score_res,score_sus = binomial_test(N_res,k_res,N_sus,k_sus)
      #plot
      PL.figure(figsize=[8,10])
      PL.title('raw scores per SNP')
      PL.subplot(311)
      PL.title('score res')
      PL.plot(pos_bin,score_res,'k.')
      PL.subplot(312)
      PL.title('score sus')
      PL.plot(pos_bin,score_sus,'k.')
      PL.subplot(313)
      PL.title('score both')
      PL.plot(pos[:,1],score,'k.')
      PL.savefig('raw_scores.png')
      #carry out simple binnin
      bin_scores,pos_bin = bin_data(score,pos,w=10E4)
      PL.figure(figsize=[8,4])
      PL.title('Binned score')
      PL.plot(pos_bin,bin_scores,'k.')
      PL.savefig('binned_scores.png')
      PL.xlim([1.9E7,2.3E7])
      PL.savefig('binned_scores_zoom.png')
    

