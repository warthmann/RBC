import scipy as SP
import pylab as PL
import os
#import cPickle
import scipy.misc
import sys
import pdb
#import pyplot

#from models import *


if __name__ == '__main__':
    X = SP.linspace(-5,5,100)

    y_res = 0.5 + 0.5*SP.exp(-SP.absolute(X))
    y_sus = 0.5 - 0.2*SP.exp(-SP.absolute(X))
    PL.figure(figsize=[10,3])
    p0=PL.plot(X,y_res,linewidth=2)
    p1=PL.plot(X,y_sus,'k-',linewidth=2)
    PL.ylim([0,1])
    PL.yticks([0,0.25,0.5,0.75,1.0])
    PL.legend([p0,p1],['res','sus'])
    PL.xlabel('Distance')
    PL.ylabel('Frequency, resistant parent')
    #PL.savefig('freq_parent.pdf')
    PL.show()
