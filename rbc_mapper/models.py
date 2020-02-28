"""module stuff
Note: the error rate for sus is currenlty completely switched off (getBeta).
TODO: think about how to get error rates that make sense.
"""

import scipy as SP
import pylab as PL
import os
import pickle
import scipy.misc
import sys
import re
import pdb
import sys
from scipy.special import gammaln, betaln
from scipy import stats


def calc_p_value_fair(data,index_res=0,p=0.5,field='counts_both'):
  """calculate the p-value of observing the res/sus ration for
  a particular observation"""
  d = data[field]
  RV = SP.zeros(d.shape[0])
  for i in xrange(d.shape[0]):
    RV[i]=stats.binom_test(d[i,index_res],d[i].sum(),p=p)
  return RV

####1. core helper functions used by all methods
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


def log_beta_binomial(k,N,alpha,beta,norm=False):
    """log pdf of the beta function
    norm: normalize? not needed if realtive ratios are considered (numerically unstable)
    """
    lpdf = 0
    if norm:
        lpdf = lnchoose(N,k)
    nom = betaln(k+alpha,N-k+beta)
    denom = betaln(alpha,beta)
    lpdf += nom
    lpdf -= denom
    return lpdf


def abtomv(X):
  """convert a,b from beta to mean and variance"""
  a = X[0]
  b = X[1]
  mean = a/(a+b)
  var  = a*b/((a+b)**2*(a+b+1))
  return [mean,var]

class RcombinationMapping(object):
  """Main class for recombination mapping"""
  def __init__(self):
      self.res =None
      self.sus =None

      self.pool_size_res = None
      self.pool_size_sus = None
      self.recombinationRate = None
      self.eps = 1E-2
      self.pos = None
      pass

  def setCountData(self,counts_res,counts_sus,pos):
      """set count data
      counts_res: [N x 2] 
      counts_sus: [N x 2]
      pos: [N x 2] (chrom, bp)
      First column needs to be resistance allele and data is required
      to be phased.
      """
      #number of res. alleles
      self.res = SP.zeros(counts_res.shape,dtype='float')
      self.sus = SP.zeros(counts_sus.shape,dtype='float')
      self.res[:,0] = counts_res[:,0]
      self.sus[:,0] = counts_sus[:,0]
      #total number of alleles
      self.res[:,1] = counts_res.sum(axis=1)
      self.sus[:,1] = counts_sus.sum(axis=1)
      #store genomic positions
      self.pos = pos
      pass

  def setPoolSizes(self,pool_size_res,pool_size_sus):
      """set Pool sizes in number of samples"""
      self.pool_size_res = pool_size_res
      self.pool_size_sus = pool_size_sus

  def setPhenotypingNoise(self,eps):
      """set rate of missphenotyping"""
      self.eps = eps
      
      
  def setRecombination(self,r):
      """recombination rate in events per bp"""
      self.recombinationRate = r

  def score(self,start_pos=None,stop_pos=None,step_size=100E3,window_size=None,opt_recombination=False,opt_eps=False):
      """stepwise scoring function
      start_pos: start position for sliding window
      stoppos  : stop  position for sliding window
      step_size: step size
      window_size: analysis window size (None). If not set, all genome-wide SNPs are jointly analzed
      opt_recombination: optimize recombination rate
      opt_eps   : optimize missphenotyping rate
      """
      if start_pos is None:
        start_pos = self.pos.min()
      if stop_pos is None:
        stop_pos = self.pos.max()
      #2. get background likelihood assuming 50:50 
      LL0 = self._LL0(SP.arange(self.res.shape[0]))

      p = start_pos
      S= []
      Sres = []
      Ssus = []
      S0 = []
      P = []
      while (p<stop_pos):
          if 1:
            #position based windowing
            dd = SP.absolute(p-self.pos)
            I  = SP.nonzero(dd<window_size)[0]
            NI = I.shape[0]
            NI = 1
          #NI = 1
          if 0:
            #total number based window
            dd = SP.absolute(p-self.pos).argmin()
            I = SP.arange(max(0,dd-100),min(self.pos.shape[0]-1,dd+100))
            NI = 1
          if 0:
            #subsample equal number of snps left and right of peak
            dd = p-self.pos
            Iw = SP.absolute(dd)<window_size
            Ip  = SP.nonzero(Iw & (dd>0))[0]
            In  = SP.nonzero(Iw & (dd<0))[0]
            Ns  = min(len(Ip),len(In),20000)
            #sample
            Irp = SP.random.permutation(len(Ip))
            Irn = SP.random.permutation(len(In))
            I = SP.concatenate((Ip[Irp][0:Ns],In[Irn][0:Ns]))
            NI = 1
            ## while True:
            ##   pdb.set_trace()
            ##   d  = self.pos[I].copy()
            ##   d[1::]-= d[0:-1]
            ##   imin = d.argmin()
            ##   if d[imin]<50:
            ##     I = SP.setdiff1d(I,I[d.argmin()])
            ##   else:
            ##     break
            ##   pass
          [score,LL_res,LL_sus]=self._LL(p,I,eps=self.eps)            
          score0 = LL0[1][I].sum() + LL0[2][I].sum()

          if opt_eps | opt_recombination:
            [score0,score] = self._LLopt(p,I,eps=self.eps,opt_eps=opt_eps,opt_recombination=opt_recombination)
         
          S.append(score/NI)
          Sres.append(LL_res.sum()/NI)
          Ssus.append(LL_sus.sum()/NI)
          S0.append(score0/NI)
          P.append(p)
          if 0:
            params_res  = self._getBeta(d=dd,pool='res')
            LL_res = self._countLL(self.res,params_res)
            mv_res      = abtomv(params_res)
            PL.ion()
            PL.figure(2,figsize=[15,6])
            PL.clf()
            print (LL_res.sum())
            PL.plot(self.pos,self.res[:,0]/self.res[:,1],'b.')
            PL.plot(self.pos,mv_res[0],'b-')
            PL.plot(self.pos,mv_res[0]+SP.sqrt(mv_res[1]),'b--')
            PL.plot(self.pos,mv_res[0]-SP.sqrt(mv_res[1]),'b--')
            pass
            
          if 0:
            PL.ion()
            PL.figure(1,figsize=[15,6])
            PL.clf()
            PL.subplot(311)
            PL.plot(self.res[:,0]/SP.double(self.res[:,1]),'b-')
            PL.plot(self.sus[:,0]/SP.double(self.sus[:,1]),'g-')
            PL.subplot(312)
            #y_res = LL_res-LL0[1][I]
            #y_sus = LL_sus-LL0[2][I]
            y_res = LL_res
            y_sus = LL_sus
            PL.plot(self.pos[I],y_res,'b-')
            PL.plot(self.pos[I],y_sus,'g-')
            PL.plot(p,0,'r*',markersize=10)
            PL.subplot(313)
            PL.plot(P,SP.array(S)-SP.array(S0),'k-')
            PL.plot(P,SP.array(Sres),'b-')
            PL.plot(P,SP.array(Ssus),'g-')
            PL.xlim([start_pos,stop_pos])
            PL.show()
            pdb.set_trace()
            pass


          #move on
          p+=step_size
      S = SP.array(S)
      S0= SP.array(S0)
      P = SP.array(P)
      return [P,S,S0]
  ### internal functions

  def _LLopt(self,p,I,ngrid=10,eps=1E-10,opt_eps=True,opt_recombination=True):
    """evluate optimized marginal likelihood:
    opt_eps: adjust rate of missphenotyping
    opt_recombination: adjust recombination rate
    """
    if opt_eps:
      epsr = 10**(SP.linspace(-10,-1,ngrid))
    else:
      epsr = [1E-2]

    if opt_recombination:
      scalemr = SP.linspace(0.1,10,ngrid)
      scalepr = SP.linspace(0.1,10,ngrid)
    else:
      scalemr = [1.0]
      scalepr = [1.0]
    LL  = SP.zeros([len(epsr),len(scalemr),len(scalepr)])
    LL0 = SP.zeros([len(epsr),len(scalemr),len(scalepr)])
    for i1 in xrange(len(epsr)):
      LL0[i1,:,:] = self._LL0(I,eps=epsr[i1])[0]
      for i2 in xrange(len(scalemr)):
        for i3 in xrange(len(scalepr)):
            LL[i1,i2,i3] = self._LL(p,I,eps=epsr[i1],scalem=scalemr[i2],scalep=scalepr[i3])[0]
    #get best parameters and return
    return [LL0.max(),LL.max()]
        
          
        
    

  def _LL(self,p,I,**kw_args):
      """evaluate likelihood score under hypothesis p is causal locus"""
      d = SP.absolute(p-self.pos[I])
      if 1:
        params_res = self._getBeta(d=d, pool='res',**kw_args)
        params_sus = self._getBeta(d=d, pool='sus',**kw_args)
      else:
        params_res = self._getBetaOLD(d=d, pool='res',**kw_args)
        params_sus = self._getBetaOLD(d=d, pool='sus',**kw_args)
      LL_res = self._countLL(self.res[I],params_res)
      LL_sus = self._countLL(self.sus[I],params_sus)
      return [(LL_res+LL_sus).sum(),LL_res,LL_sus]

  def _LL0(self,I,**kw_args):
      """evaluate likelihood score under hypothesis i  is causal locus"""
      if 1:
        params_res = self._getBeta(pool='res',podd=0.5,**kw_args)
        params_sus = self._getBeta(pool='sus',podd=0.5,**kw_args)
      else:
        params_res = self._getBetaOLD(pool='res',podd=0.5,**kw_args)
        params_sus = self._getBetaOLD(pool='sus',podd=0.5,**kw_args)
      LL_res = self._countLL(self.res[I],params_res)
      LL_sus = self._countLL(self.sus[I],params_sus)
      return [(LL_res+LL_sus).sum(),LL_res,LL_sus]


  def _podd(self,d,scalep=1.0,scalem=1.0):
      """calculate the probability of an odd number of reocmbinations
      d: distance in bp
      scalemp: scaling of recombination rate for postive and negative distances
      """
      RV = SP.zeros(d.shape)
      Ip = (d>=0)
      Im = (d<0)
      RV[Ip] = 0.5*(1-SP.exp(-2.0*d[Ip]*scalep*self.recombinationRate))
      RV[Im] = 0.5*(1-SP.exp(-2.0*d[Im]*scalem*self.recombinationRate))
      return RV

  def _Vpodd(self,d,N):
    """calculate varaince of probability of off number of recombinations
    d: distance
    N: size of pool
    """
    podd = self._podd(d)
    return 1.0/N * (1-podd)*podd
  
    
    

  def _countLL(self,counts,params,norm=False):
      """count likelihood"""
      return log_beta_binomial(counts[:,0],counts[:,1],params[0],params[1],norm=norm)

  def _meanVarBeta(self,beta_params):
      mean = beta_params[0]/(beta_params[0]+beta_params[1])
      var  = (beta_params[0]*beta_params[1]) / ( (beta_params[0]+beta_params[1])**2 (beta_params[0]+beta_params[1] +1) )
      return [mean,var]

  def _getBeta(self,d=None,podd=None,pool='res',eps=1E-2,**kw_args):
      """get parameters of beta function for distance d
      pool: res/sus
      eps: noise level from miss phenotyping
      scalep: recombination rate scaling (plus)
      scalem:  recombination rate scaling (minus)
      """
      
      #1. get podd
      if podd is None:
          podd = self._podd(d,**kw_args)
          
      if pool=='res':
        mu = 1.0-podd
        v   = 1.0/self.pool_size_res * (1-podd) * podd + eps**2

      elif pool=='sus':
        mu = 1.0/100 * (1+podd) #+ eps
        v   = 1.0/100 *  1.0/self.pool_size_sus * (1-podd) * podd + eps**2
      mu2 = mu**2
      alpha = - mu/v * (v+ mu2 - mu)
      beta  = (v + mu2 - mu)*(mu-1) / v
      return [alpha,beta]


  def _getBetaOLD(self,d=None,podd=None,pool='res',eps=1E-2):
      """get parameters of beta function for distance d
      pool: res/sus
      """
      #1. get podd
      if podd is None:
          podd = self._podd(d)
      #2. cacluclate distribution parameters
      if pool=='res':
          #factor 2 because pool size in samples, two chromosome copies
          alpha = (1-(podd+eps))*2*self.pool_size_res
          beta  = (podd+eps) * 2 * self.pool_size_res
      elif pool=='sus':
          alpha = 1.0/100 * ( 1+ (podd+eps))*2*self.pool_size_res
          beta  = 1.0/100 * (2-(podd+eps))*2*self.pool_size_sus
      #return alpha/beta
      return [alpha,beta]  

    

######old ...######

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
            PL.figure(1),PL.clf(),PL.plot(P_sus),PL.plot(P_res)
            PL.figure(2),PL.clf(),PL.plot(S)
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

    #background models
    lpdf_res0 = log_binomial(k_res,N_res,0.5)
    lpdf_sus0 = log_binomial(k_sus,N_sus,0.5)

    score_res = lpdf_res-lpdf_res0
    score_sus = lpdf_sus - lpdf_sus0
    score = score_res + score_sus

    return score,score_res,score_sus

def eval_binomial(data,index_res=0,fields=['counts_res','counts_sus'],p=[0.5,0.5]):
    #binomail test on position-wide scalee
    L = 0
    for i in xrange(len(fields)):
        field = fields[i]
        L += log_binomial(data[field][:,index_res],data[field].sum(axis=1),p[i])
    return L




def binomial_opt(data,p_res=0.99,p_sus=0.33):
    """binomail test where we pointwise optimize over the phasing information"""
    #likelihood ratios with two index_res options
    lpdf1  = eval_binomial(data,index_res=0,p=[0.99,0.33])
    lpdf2  = eval_binomial(data,index_res=1,p=[0.99,0.33])
    #background model does not requite optimization as sybmmetric
    lpdfb = eval_binomial(data,index_res=0,p=[0.5,0.5])
    lpdf  = SP.concatenate((lpdf1[:,SP.newaxis],lpdf2[:,SP.newaxis]),axis=1)
    #L
    L = lpdf.max(axis=1)
    R = {'L':L,'Lb':lpdfb,'lpdf1':lpdf1,'lpdf2':lpdf2}
    return R



             

def parse_options(argv):
    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup
    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    optional = OptionGroup(parser, 'OPTIONAL')
    required.add_option('-f', '--file', dest='data_file', metavar='FILE', help='data file with input counts', default='-')
    required.add_option('-o', '--outdir', dest='out_dir', metavar='FILE', help='output directory', default='-')
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    required.add_option('--chromosome', dest='chrom', metavar='STR', help='chromosome number', default=None)
    required.add_option('--start', dest='start', type='int', help='start bp', default=None)
    required.add_option('--stop', dest='stop', type='int', help='stop bp', default=None)
    required.add_option('--n_res', dest='n_res', type='int', help='number of samples in res pool',default=None)
    required.add_option('--n_sus', dest='n_sus', type='int', help='number of samples in sus pool',default=None)
    optional.add_option('--with_parents',  dest='with_parents', action='store_true', help='parents included in sequencing data', default=False)
    optional.add_option('--plots', dest='plots', action='store_true', help='store genome-wide figures', default=True)
    #filter options
    optional.add_option('--min-qual', dest='min_qual', type='int', help='minimal required SNP quality', default=100)
    optional.add_option('--filter-flags', dest='filter_flags', help='Comma separated list of acceptable filter flags', default='PASS')
    optional.add_option('--trust-gaTK', dest='trust_gatk', action='store_true', help='Override non-zero counts for homoczygous gatk calls', default=False)
    optional.add_option('--dp-max-ratio', dest='dp_max_ratio', type='float', help='Maximal allowed dp deviation', default=None)
    optional.add_option('--max-hs', dest='hs_max', type='float', help='Maximal allowed happlotype score', default=None)
    optional.add_option('--resolution',dest='resolution',type='int',help='Resolution in bp of the approach', default=1E5)
    optional.add_option('--min_segr_pv',dest='min_segr_pv',type='float',help='Minimal allowed segregation p-value', default=None)
    optional.add_option('--chrom_size',dest='chrom_size',type='float',help='Chromsome Size if running on a single chromosome', default=None)
    optional.add_option('--window_size',dest='window_size',type='float',help='Analysis window', default=5000E3)
    optional.add_option('--phenoNoise',dest='phenoNoise',type='float',help='Phenotyping noise', default=1E-4)
    
    #optimization flags
    
    optional.add_option('--opt_eps', dest='opt_eps', action='store_true', help='verbosity', default=False)
    optional.add_option('--opt_recombination', dest='opt_recombination', action='store_true', help='verbosity', default=False)

    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) < 3:
        parser.print_help()
        sys.exit(2)

    return options


   
  

def trivial_phasing(data,max_coverage=0,do_phase=True):
    """preprocess and phase the data:
    major allele in res is k"""

    counts_res = data['counts_res']
    counts_sus = data['counts_sus']
    #separate in N and k
    #k is the major allele in res!
    iN = SP.arange(counts_res.shape[0])
    if do_phase:
        ik = counts_res.argmax(axis=1)
    else:
        ik = SP.zeros(counts_res.shape[0],dtype='int')
    #get total count and k for res and sus
    N_res = counts_res.sum(axis=1)
    k_res = counts_res[iN,ik]
    N_sus = counts_sus.sum(axis=1)
    k_sus = counts_sus[iN,ik]
    return N_res,k_res,N_sus,k_sus
