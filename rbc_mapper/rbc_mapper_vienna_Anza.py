import scipy as SP
import pylab as PL
import os
import pickle
import scipy.misc
import sys
import pdb
import copy 
import yaml
import numpy as np
import logging

from models_Anza import *
from data_Anza import *



#from old.models import *
#import old.models as OM


def analyse_rm(data,samples_res,samples_sus,recombination_rate,window_size,step_size,opt_eps=False,opt_recombination=False,phenoNoise=0,out_base=None,EMAR_SusPool=0):
    rm = RcombinationMapping()
    rm.setCountData(data['counts_res'],data['counts_sus'],data['pos'],EMAR_SusPool)
    rm.setPoolSizes(samples_res,samples_sus)
    rm.setRecombination(recombination_rate)
    rm.setPhenotypingNoise(phenoNoise)
    [P,S,S0]=rm.score(step_size=step_size,window_size=window_size,opt_eps=opt_eps,opt_recombination=opt_recombination)

    score = S-S0
    #set scores to 0 if negative
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
        PL.grid()
        PL.title("Chromosome : %s" %np.unique(data['chrom']))
        PL.show()

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
        M= np.zeros([P.shape[0],3],dtype='object')
        M[:,0] = data['chrom'][0]
        M[:,1] = P
        M[:,2] = score
        np.savetxt(out_base+'.csv',M,fmt='%s',delimiter=',')
        



def main(options):
    
    #2. hard coded preprocessing parameters
    preprocess_params = {'min_qual':options["min_qual"],'trust_gatk':options["trust_gatk"],'dp_max_ratio':options["dp_max_ratio"],'hs_max':options["hs_max"],'filter_flags':options["filter_flags"],'chrom':options["chrom"],'start':options["start"],'stop':options["stop"],'min_segr_pv':options["min_segr_pv"]}
    analyse_params = {'step_size':options["resolution"],'window_size': options["window_size"],'opt_eps':options["opt_eps"],'opt_recombination':options["opt_recombination"],'phenoNoise':options["phenoNoise"],'EMAR_SusPool': options["EMAR_SP"]}

    samples_res = options["n_res"]
    samples_sus = options["n_sus"]
    data_file = options["data_file"]
    out_dir  = options["out_dir"]
    #data file and caching file:
    data_file_cached = data_file +'.pickle'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if samples_res is None:
        print ("Need number of res samples")
        sys.exit(1)
    if samples_sus is None:
        print ("Need number of sus samples")
        sys.exit(1)
# if recalc parameter is given on the command line a new pickle file will be generated
    recalc = 'recalc' in sys.argv # "recalc" in commandline triggers replacement of existing pickle file 
    #1. load raw data

# The below code is to get the data from the data_file or if there is a pickle file use that.
    fileName, fileExtension = os.path.splitext(data_file)
    if not os.path.exists(data_file_cached) or recalc:
        data = load_data_irbc(data_file,with_parents=options["with_parents"])
        pickle.dump(data,open(data_file_cached,'wb'),-1)
    else:
        data = pickle.load(open(data_file_cached,'rb'))

# The below code is for checking if the chromosome names are defined in the command line or config file if not then it will ask the user to either enter the required chromosomes or run for all the chromosomes(Anza)

    if preprocess_params['chrom'] is None:
        print('\n\n\nWarning: All the chromosomes will be used:%s' %np.unique(data['chrom']))
        print('\nWarning: Total number of chromosomes: %s' %len(np.unique(data['chrom'])))
        chromosomes_selected= input('\n\nPress Enter to use all the chromosomes or the Names that you want comma seprated:\n')
        if chromosomes_selected is None:
            chromosome_names= np.unique(data['chrom'])
        else:
            chromosome_names= str.split(chromosomes_selected,',')
    else:
        chromosome_names = str.split(preprocess_params['chrom'],',')
# Once the chromosome names are defined or given the below code will send chromosome name one by one for processing(This is done so that the original code is changed).(Anza)
    #filter_flags_split= str.split(preprocess_params['filter_flags'],'@')
    for i in range(0,len(chromosome_names)):
        print('\n\n\nChromosome: %s' %chromosome_names[i])
        print('Filter Flags: %s' %preprocess_params['filter_flags'])

    #1. load raw data
        fileName, fileExtension = os.path.splitext(data_file)
        if not os.path.exists(data_file_cached) or recalc:
            data = load_data_irbc(data_file,with_parents=options["with_parents"])
            pickle.dump(data,open(data_file_cached,'wb'),-1)
        else:
            data = pickle.load(open(data_file_cached,'rb'))

        chrom=chromosome_names[i]
        preprocess_params['chrom']=chromosome_names[i]
        #preprocess_params['filter_flags']=filter_flags_split[i]
# The statement is for sending the data and the parameters to the preprocess_data function where all the filtering is done wrt the parameters set and we get back the filtered file.(Anza)
        preprocess_data(data,**preprocess_params)    
# The below statement is to get a parameter which should be multiplied with the counts_res and counts_sus to make the weight of these two approximately same for counts both calculation.(Anza)      
        #get library size correction
        [LSres,LSsus] = lib_size_factors(data)
        

        ## Anza Testing
        ##[LSres,LSsus] = lib_size_factors(data)


# counts_both is adding allel count for res and sus. Right now the counts_res has major minor configuration and counts_sus has a configuration in which major or minor allel can be anywhere.
        data['counts_both'] = LSres*data['counts_res']+LSsus*data['counts_sus']
    #    data['counts_both'] = LSres*data['counts_res']+LSsus*data['counts_sus']
        #print(LSres*data['counts_res']) #printed for better understanding (Norman, Anza Jan 2020)
        #print(LSsus*data['counts_sus'])
        #print(data['counts_both'])
        
        #filter on allele counts 
        if 1:
            #filter by ratio (a normal SNP should have 50% over both pools, but this is true only, if the numbers of individuals in th epools reflect the inheritance (I.e. 3:1, in a recessive one-gene trait, the phenotype(-) pools has to have 3 times the size of the phenotype(+) pool)
            # this filter needs to take pools size into account. Do not hard code ratios.
            ratio =data['counts_both'][:,0]/data['counts_both'].sum(axis=1)
            #Iok = (0.45<ratio) & (ratio<0.65)
            Iok = (0.30<ratio) & (ratio<0.70)
            filter_data(data,Iok)

        pos = data['pos']
        chrom = data['chrom']

        #loop over chromosomes and carry out analysis
        uchrom = np.unique(chrom)

        for c in uchrom:
            
            _data = copy.deepcopy(data)
            Ic = (data['chrom']==c)
            filter_data(_data,Ic)
            
            # changed none to 0 in the below expression by Anza
            if options["chrom_size"] is None:
                recombination_rate = options["EX_CO_Chrom"]/(data['pos'].max()-data['pos'].min())
                
            else:
                
                recombination_rate = options["EX_CO_Chrom"]/options["chrom_size"]

  # The below if statement is not working 
            if 0:
                print ("recombination rate factor on")
                recombination_rate *= 2

    # This is where the plot is being developed
            analyse_rm(_data,samples_res,samples_sus,recombination_rate,out_base = os.path.join(out_dir,'chrom_%s' % (c)),**analyse_params)


        if 1:
            #PL.figure()
            PL.subplot(411)
            PL.title("counts_both")
            PL.plot(data['pos'],data['counts_both'][:,0]/data['counts_both'].sum(axis=1),'b.',alpha=0.2)
            PL.ylim(-0.1, 1.1)
            PL.grid()
            PL.subplot(412)
            PL.title("counts determined major allele in res pool")
            PL.plot(data['pos'],np.array(data['counts_res'][:,0],dtype='float')/data['counts_res'].sum(axis=1),'b.',alpha=0.2)
            PL.ylim(-0.1, 1.1)
            PL.grid()
            PL.subplot(413)
            PL.title("counts this same allele in sus pool")
            PL.plot(data['pos'],np.array(data['counts_sus'][:,0],dtype='float')/data['counts_sus'].sum(axis=1),'b.',alpha=0.2)
            PL.ylim(-0.1, 1.1)
            PL.grid()
            PL.subplot(414)
            PL.title("counts determined minor allele in res pool")
            PL.plot(data['pos'],np.array(data['counts_res'][:,1],dtype='float')/data['counts_res'].sum(axis=1),'b.',alpha=0.2)
            PL.ylim(-0.1, 1.1)
            PL.grid()
            PL.show()


if __name__ == '__main__':
    #parse command line
    file_name= parse_options(sys.argv)
    
    logging.basicConfig(filename="Information.log",format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO,datefmt='%Y-%m-%d %H:%M:%S')  

# The below statements are for checking if the parameters are set on the command line or if there is a config file from where the data should be taken(Anza)
    if len(sys.argv) > 3:
        options_parsed   = parse_options(sys.argv)
        options = vars(options_parsed)
        print(options)
        logging.info('The parameters are set as: %s' %options)
        Command_line_argv_file = 'argv.txt'
        f= open(Command_line_argv_file,'w')
        f.write(str(options).replace(',','\n'))
    else:
        options_parsed= vars(parse_options(sys.argv))
        file_name=options_parsed['config_file']
        print('Config_file being used: %s' %file_name)
        
   
        with open(file_name) as file:
            options = yaml.load(file)
            print(options)
        
    
        
    #print (type(options))
    #print (options['min_qual'])
    '''with open('config_1.yaml') as file:
        options = yaml.load(file)
        print(options)'''

    #print ("Note: using options: %s" % (str(options)))
    #print(type(options['window_size']))
    #print(type(options['window_size']))
    logging.info('The parameters are set as: %s' %options)
    main(options)
    

