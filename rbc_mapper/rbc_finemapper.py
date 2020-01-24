import scipy as SP
import pylab as PL
import os
import cPickle
import scipy.misc
import sys
import pdb
import copy 

from models import *
from data import *


if __name__ == '__main__':
    #parse command line
    options   = parse_options(sys.argv)
    print "Note: using options: %s" % (str(options))
    rel_score = True

    xc = 20129204
    
    #2. hard coded preprocessing parameters
    preprocess_params = {'min_qual': options.min_qual,'trust_gatk':options.trust_gatk,'dp_max_ratio':options.dp_max_ratio,'hs_max':options.hs_max,'filter_flags':options.filter_flags,'chrom':options.chrom,'start':options.start,'stop':options.stop,'min_segr_pv':options.min_segr_pv}
    analyse_params = {'step_size':options.resolution,'window_size': options.window_size,'opt_eps':options.opt_eps,'opt_recombination':options.opt_recombination,'phenoNoise':options.phenoNoise}
        
    samples_res = options.n_res
    samples_sus = options.n_sus
    data_file = options.data_file
    out_dir  = options.out_dir
    #data file and caching file:
    data_file_cached = data_file +'.pickle'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if samples_res is None:
        print "Need number of res samples"
        sys.exit(1)
    if samples_sus is None:
        print "Need number of sus samples"
        sys.exit(1)
   
    recalc = 'recalc' in sys.argv
    #1. load raw data
    fileName, fileExtension = os.path.splitext(data_file)
    if not os.path.exists(data_file_cached) or recalc:
        data = load_data_irbc(data_file,with_parents=options.with_parents)
        cPickle.dump(data,open(data_file_cached,'wb'),-1)
    else:
        data = cPickle.load(open(data_file_cached,'rb'))

    preprocess_data(data,**preprocess_params)    

    #get library size correction
    [LSres,LSsus] = lib_size_factors(data)
    data['counts_both'] = LSres*data['counts_res']+LSsus*data['counts_sus']

    if 0:
        #filter ratio
        ratio =data['counts_both'][:,0]/data['counts_both'].sum(axis=1)
        Iok = (0.45<ratio) & (ratio<0.8)
        filter_data(data,Iok)

    #chromosome length:
    pos = data['pos']
    chrom = data['chrom']
    globals().update(data)

    cf = (counts_both[:,0]/counts_both.sum(axis=1))/0.66

    #1. initialize recombination mapping
    recombination_rate = 1.0/(2.71E7)
    rm = RcombinationMapping()
    rm.setCountData(data['counts_res'],data['counts_sus'],data['pos'])
    rm.setPoolSizes(samples_res,samples_sus)
    rm.setRecombination(recombination_rate)
    rm.setPhenotypingNoise(analyse_params['phenoNoise'])
   
    if 1:
        #using fitter and apply to these data
        [P,S,S0]=rm.score(step_size=100E3,window_size=5E7,opt_eps=analyse_params['opt_eps'],opt_recombination=analyse_params['opt_recombination'])
        PL.plot(P,S,'k-')
        PL.savefig(os.path.join(out_dir,'mapping.pdf'))

    if 1:
        PL.figure()
        #0. plot theoretical curve around xc
        pos_range = SP.linspace(pos.min(),pos.max(),1000)
        D_range   = SP.absolute(pos_range-(xc+0.01E7))
        podd     = rm._podd(D_range)
        Spodd    = SP.sqrt(rm._Vpodd(D_range,options.n_res))
        #1. plot theory
        rt = 0.98
        PL.plot(pos_range,rt-podd,'k-')
        PL.plot(pos_range,(rt-podd) + Spodd,'k--')
        PL.plot(pos_range,(rt-podd) - Spodd,'k--')
        #1. plot raw data
        PL.plot(pos,SP.double(counts_res[:,0])/counts_res.sum(axis=1),'b.')
        PL.savefig(os.path.join(out_dir,'fit.pdf'))
        
    if 0:

        PL.figure()
        PL.subplot(311)
        PL.plot(pos,SP.double(counts_sus[:,0])/counts_sus.sum(axis=1),'b.')
        PL.subplot(312)
        PL.plot(pos,SP.double(counts_res[:,0])/counts_res.sum(axis=1),'b.')
        PL.subplot(313)
        PL.plot(pos,SP.double(counts_both[:,0])/counts_both.sum(axis=1),'r.')
    
    
