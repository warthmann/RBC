import scipy as SP
import pylab as PL
import os
import pickle
import scipy.misc
import sys
import re
import pdb
import string
import sys
from models import calc_p_value_fair



def load_data_irbc(data_file,with_parents=False):
    """load fixed format vcf file used with the scripts provided here
    with_parents: parents included in data file?
    """
    #0. define regular expressions
    m_comment = re.compile('#.*')
    #1. create result fields
    pos = []
    chrom = []
    counts_res = []
    counts_sus = []
    alleles = []
    geno_res = []
    geno_sus = []
    qual = []
    filter = []
    type = []
    segr = []
    hapsc = []
    i_p0 = None
    i_p1 = None
    parent_geno = None

    if with_parents:
        i_res = 10
        i_sus = 11
        i_p0  = 8
        i_p1 = 9
        parent_geno = []
    else:
        i_res = 13
        i_sus = 14
        
    #2. open file
    f = open(data_file,'rU')
    #3. parse
    #3.1 skip comment lines
    while True:
        line=f.readline()
        if not m_comment.match(line):
            break
    while True:
        #read next line
        line = f.readline()
        ls = str.split(line,sep='\t')
        #split successful?
        if not ls:
            break
        if len(ls)<2:
            break
        #segr
        _segr = ls[10]
        if _segr=='NA':
            continue
        #type filter
        _type = ls[9]
        if _type!='SNP':
            continue
        #parse line
        #positional information
        _chrom = ls[0]
        _pos   = int(ls[1])
        chrom.append(_chrom)
        pos.append(_pos)
        #alleles
        _alleles = ls[3:7]
        alleles.append(_alleles)
        #qual
        qual.append(float(ls[7]))
        #filter flag
        #_filter = string.upper(ls[8])
        _filter = ls[8].upper()
        #split
        _filter = str.split(_filter,',')
        filter.append(_filter)
        type.append(_type)
        segr.append(_segr)
        _hapsc = ls[11]
        #hapsc_ = string.replace(_hapsc,'NA','NAN')
        hapsc_ = _hapsc.replace('NA','NAN')
        hapsc.append(hapsc_)
        #parent information
        if i_p0 is not None:
            p0 = str.split(ls[i_p0],':')[0]
            p1 = str.split(ls[i_p1],':')[0]
            parent_geno.append([p0,p1])
        #3. counts
        ## description string
        descr = ls[12]
        res = ls[i_res]
        sus = ls[i_sus]
        descr_fields = str.split(descr,':')
        res_fields   = str.split(res,':')
        sus_fields   = str.split(sus,':')
        for i in range(len(descr_fields)):
            field = descr_fields[i]
            if field=='GT':
                geno_res.append(res_fields[i])
                geno_sus.append(sus_fields[i])
            elif field=='AD':
                counts_res.append(SP.array(str.split(res_fields[1],','),dtype='uint16'))
                counts_sus.append(SP.array(str.split(sus_fields[1],','),dtype='uint16'))
        pass
    if parent_geno is not None:
        parent_geno = SP.array(parent_geno)
    RV = {'pos':SP.array(pos),'chrom':SP.array(chrom), 'counts_res':SP.array(counts_res),'counts_sus':SP.array(counts_sus),'alleles_res':SP.array(alleles),'alleles_sus': SP.array(alleles), 'qual':SP.array(qual,dtype='int'),'parent_geno':parent_geno,'qual':SP.array(qual,dtype='float'),'type':SP.array(type),'segr':SP.array(segr),'hapsc':SP.array(hapsc,dtype='float'),'filter':SP.array(filter,dtype='object'),'geno_res':SP.array(geno_res),'geno_sus':SP.array(geno_sus)}
    return RV


def load_data_vcf(data_file,with_parents=False):
    """load fixed format vcf file used with the scripts provided here
    with_parents: parents included in data file?
    """
    #0. define regular expressions
    m_comment = re.compile('#.*')
    m_chrom   = re.compile('.*chr(.*)')
    #1. create result fields
    pos = []
    chrom = []
    counts_res = []
    counts_sus = []
    alleles = []
    geno_res = []
    geno_sus = []
    qual = []
    hscore = []
    flag = []
    i_p0 = None
    i_p1 = None
    parent_geno = None

    if with_parents:
        i_res = 10
        i_sus = 11
        i_p0  = 8
        i_p1 = 9
        parent_geno = []
    else:
        i_res = 8
        i_sus = 9
        

    #2. open file
    f = open(data_file,'rU')
    #3. parse
    #3.1 skip comment lines
    while True:
        line=f.readline()
        if not m_comment.match(line):
            break
    while True:
        ls = str.split(line,sep='\t')
        #split successful?
        if not ls:
            break
        if len(ls)<2:
            break
        #parse line
        #positional information
        _chrom = m_chrom.match(ls[0]).group(1)
        _pos   = int(ls[1])
        chrom.append(_chrom)
        pos.append(_pos)
        #alleles
        _alleles = ls[3:5]
        #qual
        qual.append(float(ls[5]))
        alleles.append(_alleles)
        #flat
        flag.append(ls[6])
        #parent information
        if i_p0 is not None:
            p0 = str.split(ls[i_p0],':')[0]
            p1 = str.split(ls[i_p1],':')[0]
            parent_geno.append([p0,p1])
        #3. counts
        ## description string
        descr = ls[7]
        res = ls[i_res]
        sus = ls[i_sus]
        descr_fields = str.split(descr,':')
        res_fields   = str.split(res,':')
        sus_fields   = str.split(sus,':')
        for i in range(len(descr_fields)):
            field = descr_fields[i]
            if field=='GT':
                geno_res.append(res_fields[i])
                geno_sus.append(sus_fields[i])
            elif field=='AD':
                counts_res.append(SP.array(str.split(res_fields[1],','),dtype='uint16'))
                counts_sus.append(SP.array(str.split(sus_fields[1],','),dtype='uint16'))
        #read next line
        line = f.readline()
        pass
    if parent_geno is not None:
        parent_geno = SP.array(parent_geno)
    RV = {'pos':SP.array(pos),'chrom':SP.array(chrom), 'counts_res':SP.array(counts_res),'counts_sus':SP.array(counts_sus),'alleles_res':SP.array(alleles),'alleles_sus': SP.array(alleles), 'qual':SP.array(qual,dtype='int'),'parent_geno':parent_geno,'filter':SP.array(filter)}
    return RV


def preprocess_data(D,filterN=True,filter_dash=True,enforce_match_major=True,enforce_match_minor=False, res_index=0,filter_parent=True,min_qual=0,trust_gatk=False,filter_flags='PASS',dp_max_ratio=0.5,hs_max=1.0,chrom=None,start=None,stop=None,min_segr_pv=None):
    """preprocess loaded data file
    D: data object from rbci reader
    filterN: filter alleles with N anywhere
    filter_dash: filter alleles with - anywhere
    enforce_match_major: enforce major allele match between pools
    enforce_match_minor: enforce minor allele match between pools
    res_index: index of resistant allelle (default: major =0)
    min_qual: minimum required quality
    trust_gatk: set minor count to zero if hom call
    filter_flag: comma separated list with legal filter flags
    dp_max_ratio: maximum relative deviation of dp max
    hs_max: max happlotype score
    chrom/start/stop: window where da is being preprocessed.
    """
    #1.calculate coverage averages across the genome
    mCres = SP.median(D['counts_res'].sum(axis=1))
    mCsus = SP.median(D['counts_sus'].sum(axis=1))

    #2.filter genomic positiosn if applicable
    Iok = SP.ones(D['alleles_res'].shape[0],dtype='bool')
    if chrom is not None:
        Iok = Iok & (D['chrom']==chrom)
    if start is not None:
        Iok = Iok & (D['pos']>=start)
    if stop is not None:
        Iok = Iok & (D['pos']<=stop)

    print ("Note: restricting analysis to chrom: %s, start: %s, stop: %s" % (chrom,start,stop))

    #filter
    filter_data(D,Iok)
    
    #do all other filters
    alleles_res = D['alleles_res']
    alleles_sus = D['alleles_sus']
    counts_res  = D['counts_res']
    counts_sus  = D['counts_sus']
    Iok = SP.ones(alleles_res.shape[0],dtype='bool')    

    if filterN:
      #1. remove consistent "N" in either res or sus:
      Iok  = Iok  & ((D['alleles_res']!='N').any(axis=1) & (D['alleles_sus']!='N').any(axis=1))
    if filter_dash:
      #1. remove any '-' in res or sus:
      Iok  = Iok  & (D['alleles_res']!='-').all(axis=1) & (D['alleles_sus']!='-').all(axis=1)
    #enforce match major?
    if enforce_match_major:
      Iok = (alleles_res[:,0]==alleles_sus[:,0])
    #enforce match minor?
    if enforce_match_minor:
      Iok = Iok & (alleles_res[:,1]==alleles_sus[:,1])
    #filter quality
    if min_qual:
        Iok = Iok & (D['qual'].squeeze()>=min_qual)
    #filter hs score
    if hs_max:
        Iok = Iok & (D['hapsc']<=hs_max)
        pass
    #filter filter flag
    #1. get list of legal flags
    if filter_flags:
        #lf = str.split(string.upper(filter_flags),',')
        lf = str.split(filter_flags.upper(),',')
        ifilter = SP.array([SP.array([e in lf for e in _filter]).all() for _filter in D['filter']])
        Iok = Iok & ifilter
    
    #coverage filter
    #res
    Cres = D['counts_res'].sum(axis=1)
    Cs = dp_max_ratio*mCres
    c_min = mCres-Cs
    c_max = mCres+Cs
    Iok = Iok & (Cres<c_max) & (Cres>c_min)
    #sus
    Csus = D['counts_sus'].sum(axis=1)
    Cs = dp_max_ratio*mCsus
    c_min = mCsus-Cs
    c_max = mCsus+Cs
    Iok = Iok & (Cres<c_max) & (Cres>c_min)
    #filter parent geno for segragation if we do have parental data
    if filter_parent and D['parent_geno'] is not None:
        ip = (D['parent_geno'][:,0]=='0|0') & (D['parent_geno'][:,1]=='1|1') | ((D['parent_geno'][:,0]=='1|1') & (D['parent_geno'][:,1]=='0|0'))
        Iok = Iok & ip
        pass
    #apply filter
    print ("Note: filtering retained %d/%d SNPs" % (Iok.sum(),Iok.shape[0]))
    filter_data(D,Iok)
    #calc joint count with lib size factor adjustment
    [LSres,LSsus] = lib_size_factors(D)
    D['counts_both'] = LSres*D['counts_res']+LSsus*D['counts_sus']
    
    #segragation p-value?
    if min_segr_pv is not None:
        pv = calc_p_value_fair(D)
        Iok = pv>min_segr_pv
        filter_data(D,Iok)
    
    #gatk transform
    if trust_gatk:
        #get positiosn where homozgyous call was made and set counts to 0 (minor)
        Ibadr = (D['geno_res']=='0|0') | (D['geno_res']=='1|1') | (D['geno_res']=='1/1') | (D['geno_res']=='0/0')
        D['counts_res'][Ibadr,1] = 0
        Ibads = (D['geno_sus']=='0|0') | (D['geno_sus']=='1|1') | (D['geno_sus']=='1/1') | (D['geno_sus']=='0/0')
        D['counts_sus'][Ibads,1] = 0                                      
    pass



def bin_data(d,pos,window=None,step=20E3,average=False,N=None):
    # 1. binning
    d_bin = []
    pos_bin = []
    if window is not None:
        p0 = 0
        pmax = pos.max()
        while p0 < pmax:
            p1 = p0 + window
            Ib = (p0<=pos) & (pos<p1)
            p0+= step
            if not Ib.any():
              continue
            dd = SP.array(d[Ib,:].sum(axis=0),dtype='float')
            if average:
                dd/= Ib.sum()
            d_bin.append(dd)
            pos_bin.append(pos[Ib][0])
    elif N is not None:
        i0 = 0
        while i0<d.shape[0]:
            i1 = i0 +N
            _d = d[i0:i1].sum()
            _p = pos[i0:i1].mean()
            d_bin.append(_d)
            pos_bin.append(_p)
            i0+=N
    d_bin = SP.array(d_bin)
    pos_bin    = SP.array(pos_bin)
    return d_bin,pos_bin


def filter_data(D,I):
    """subset all elements in data with index I"""
    for k in D.keys():
        if D[k] is not None:
            D[k] = D[k][I]




def lib_size_factors(data):
    """calculate library size correction factors"""
    res_sum = data['counts_res'].sum()
    sus_sum = data['counts_sus'].sum()
    print(res_sum)
    print(sus_sum)

    #L = [1.0, SP.double(res_sum/sus_sum)] #corrected the direction for normalisation (Norman and Anza, jan 2020)
    L = [SP.double(res_sum/sus_sum), 1.0]
    #print(L)
    return L
         
    
