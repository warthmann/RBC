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



#from old.models import *
#import old.models as OM


def analyse_rm(data,samples_res,samplse_sus,recombination_rate,window_size,step_size,opt_eps=False,opt_recombination=False,phenoNoise=0,out_base=None):
    rm = RcombinationMapping()
    rm.setCountData(data['counts_res'],data['counts_sus'],data['pos'])
    rm.setPoolSizes(samples_res,samples_sus)
    rm.setRecombination(recombination_rate)
    rm.setPhenotypingNoise(phenoNoise)
    [P,S,S0]=rm.score(step_size=step_size,window_size=window_size,opt_eps=opt_eps,opt_recombination=opt_recombination)

    score = S-S0
    score[score<0] = 0
    if 0:
        PL.figure(1,figsize=[10,5])
        PL.clf()
        PL.subplot(211)
        PL.plot(P,S0,'k--')
        PL.plot(P,S)
        PL.subplot(212)
        PL.plot(P,score)
    if 1:
        PL.figure(1,figsize=[10,5])
        PL.clf()
        PL.plot(P,score,'k.',alpha=0.5)
        #PL.show()
    if 0:
        PL.figure(1,figsize=[10,5])
        PL.clf()
        PL.subplot(211)
        PL.plot(P,score,'k.',alpha=0.5)
        #plot quality of segragation assumption
        PL.subplot(212)
        L0 = calc_p_value_fair(data)
        L0_bin,p = bin_data(L0,data['pos'])
        PL.plot(data['pos'],-SP.log10(L0),'b.')
        PL.ion()
        PL.show()
        pass

    #output?
    if out_base is not None:
        #1. figure
        PL.savefig(out_base+'.png')
        #2. csv file with score
        M= SP.zeros([P.shape[0],3],dtype='object')
        M[:,0] = data['chrom'][0]
        M[:,1] = P
        M[:,2] = score
        SP.savetxt(out_base+'.csv',M,fmt='%s',delimiter=',')
        



if __name__ == '__main__':
    #parse command line
    options   = parse_options(sys.argv)
    print "Note: using options: %s" % (str(options))

    
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
    #chromosome length:

    if 0:
        #filter ratio
        ratio =data['counts_both'][:,0]/data['counts_both'].sum(axis=1)
        Iok = (0.45<ratio) & (ratio<0.65)
        filter_data(data,Iok)

    pos = data['pos']
    chrom = data['chrom']

    #loop over chromosomes and carry out analysis
    uchrom = SP.unique(chrom)

    for c in uchrom:
        print "Note: processing Chromosome %s" % (c)
        _data = copy.deepcopy(data)
        Ic = (data['chrom']==c)
        filter_data(_data,Ic)
        if options.chrom_size is None:
            recombination_rate = 1.0/(data['pos'].max()-data['pos'].min())
        else:
            recombination_rate = 1.0/options.chrom_size
        if 0:
            print "recombination rate factor on"
            recombination_rate *= 1
        analyse_rm(_data,samples_res,samples_sus,recombination_rate,out_base = os.path.join(out_dir,'chrom_%s' % (c)),**analyse_params)


    if 1:
        PL.figure()
        PL.subplot(311)
        PL.plot(data['pos'],data['counts_both'][:,0]/data['counts_both'].sum(axis=1),'b.',alpha=0.2)
        PL.subplot(312)
        PL.plot(data['pos'],SP.array(data['counts_res'][:,0],dtype='float')/data['counts_res'].sum(axis=1),'b.',alpha=0.2)
        PL.subplot(313)
        PL.plot(data['pos'],SP.array(data['counts_sus'][:,0],dtype='float')/data['counts_sus'].sum(axis=1),'b.',alpha=0.2)

