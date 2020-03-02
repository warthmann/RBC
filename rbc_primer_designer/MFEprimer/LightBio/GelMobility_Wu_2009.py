#!/usr/bin/env python
# -*- coding: utf-8 -*- #   
'''Mobility of the PCR amplicons in agarose gel electrophoresis
was calculated following the semi-log formula reported by Helling
et al [Ref. 1] in 1974.

    Y = a - b * ln (X + k)

  Where Y is the mobility and X is the size of the PCR amplicons.

  To get an reliable prediction of mobility of PCR amplicons, we
determined these parameters (a, b and k) by statistically analyzing
the commercially available DNA markers of known size (100 bp, 250 bp,
500 bp, 750 bp, 1000 bp and 2000 bp, from TaKaRa Inc.) in agarose gel
with different concentrations. The DNA fragment of 100 bp was taken as
a reference band in the agarose gel electrophoresis, and the relative
mobility was determined according to the 100 bp DNA band. When the 
mobility of the reference band reached 5 cm in running the agarose gel,
the electrophoresis was stoped and the mobility of other DNA bands
were measure to determine the relative mobility. Briefly, the 
parameters in the abovementioned formula were determined on the 
0.5% ~ 2% agarose gel electrophoresis with triplicate measurements
in out lab and shown below for producing the virtual electrophoretogram
in MPprimer program.

  The parameters obtained in our lab for Helling's formula:

0.5% (Gel concentration)

    a = 2.7094
    b = 0.2691
    k = 464.4412

1.0% (Gel concentration)

    a = 2.3977
    b = 0.2700
    k = 73.9788

1.5% (Gel concentration)

    a = 2.3221
    b = 0.2634
    k = 48.0873

2.0% (Gel concentration)

    a = 2.1333
    b = 0.2561
    k = 18.5417

  The parameters obtained in our lab for Yonghong Wu's formula:

0.5% (Gel concentration)

    a = 3.481783
    b = 0.4968585
    k = 1001.76

1.0% (Gel concentration)

    Y=exp(4.783911-0.7481319*ln(X+496.7385))      

    a = 4.783911
    b = 0.7481319
    k = 496.7385

1.5% (Gel concentration)

    Y=exp(4.607372-0.7360514*ln(X+420.3904))      

    a = 4.607372
    b = 0.7360514
    k = 420.3904

2.0% (Gel concentration)

    Y=exp(5.46549-0.9068884*ln(X+311.2767))       

    a = 5.46549
    b = 0.906884
    k = 311.2767

Reference

[1] Helling, R. B., Goodman, H. M., & Boyer, H. W. (1974). Analysis of 
endonuclease R-EcoRI fragments of DNA from lambdoid bacteriophages and 
other viruses by agarose-gel electrophoresis. Journal of virology, 
14(5), 1235-44.
[2] Yonghong Wu et al. 2009, under review.
'''

# Author: Wubin Qu [quwubin@gmail.com], BIRM, China
# Date: 2009-10-12
# License: GPL v3

import sys

def load_gel_para_dict(gel_conc=1.0, formula='Wu'):
    gel_para_dict = {}
    gel_para_dict['Helling'] = {}
    gel_para_dict['Wu'] = {}

    # Parameter for Helling's formula
    gel_para = {}
    gel_para[0.5] = {
	'a' : 2.709423,
	'b' : 0.2691262,
	'k' : 464.4412,
    }

    gel_para[1.0] = {
	'a' : 2.397726,
	'b' : 0.2700602,
	'k' : 73.97887,
    }

    gel_para[1.5] = {
	'a' : 2.322184,
	'b' : 0.2634794,
	'k' : 48.08739,
    }

    gel_para[2.0] = {
	'a' : 2.133323,
	'b' : 0.2561728,
	'k' : 18.54176,
    }

    gel_para_dict['Helling'] = gel_para

    # Parameter for Wu's formula
    gel_para = {}

#0.5% (Gel concentration)
    gel_para[0.5] = {
        #Y=exp(3.481783-0.4968585*ln(X+1001.76))      
        'a' : 3.481783,
        'b' : 0.4968585,
        'k' : 1001.76,
    }

#1.0% (Gel concentration)
    gel_para[1.0] = {
        #Y=exp(4.783911-0.7481319*ln(X+496.7385))      
        'a' : 4.783911,
        'b' : 0.7481319,
        'k' : 496.7385,
    }

#1.5% (Gel concentration)
    gel_para[1.5] = {
        #Y=exp(4.607372-0.7360514*ln(X+420.3904))      
        'a' : 4.607372,
        'b' : 0.7360514,
        'k' : 420.3904,
    }

#2.0% (Gel concentration)
    gel_para[2.0] = {
        #Y=exp(5.46549-0.9068884*ln(X+311.2767))       
        'a' : 5.46549,
        'b' : 0.906884,
        'k' : 311.2767,
    }

    gel_para_dict['Wu'] = gel_para

    err_msg = 'Gel concentration is illegal, Currently, only 0.5, 1, 1.5, 2.0 are allowed.'

    try: 
        gel_conc = float(gel_conc)
        if formula == 'Wu':
            a = gel_para_dict['Wu'][gel_conc]['a']
            b = gel_para_dict['Wu'][gel_conc]['b']
            k = gel_para_dict['Wu'][gel_conc]['k']
        else:
            a = gel_para_dict['Helling'][gel_conc]['a']
            b = gel_para_dict['Helling'][gel_conc]['b']
            k = gel_para_dict['Helling'][gel_conc]['k']
    except:
	print >> sys.stderr, err_msg
	exit()

    return gel_para_dict, a, b, k

def get_size_range(size, gel_conc=1.0, ref_mobility=50, offset=2, formula='Wu'):
    Y = cal_mobility(size, gel_conc)
    # Set 2 mm as the distance which the bands can be \
            #seperated by naked eyes
    Xmin = cal_size(Y + offset, gel_conc=gel_conc, ref_mobility=ref_mobility, formula='Wu')
    Xmin = cal_size(Y - offset, gel_conc=gel_conc, ref_mobility=ref_mobility, formula='Wu')

    return Xmin, Xmax

def cal_mobility(X, gel_conc=1.0, ref_mobility=50, formula='Wu'):
    '''Cal mobility based on size'''
    import math
    gel_para_dict, a, b, k = load_gel_para_dict(gel_conc=gel_conc, formula=formula)

    X = float(X)
    gel_conc = float(gel_conc)

    # X: size (bp)
    # ref_mobility: the mobility distance of the fastest DNA segment

    if formula == 'Wu':
        Y = math.exp(a - b * math.log(X + k))
    else:
        Y = a - b * math.log(X + k)

    # Y: the relative mobility = mobility distance / ref_mobility
    Y = '%.1f' % (Y * ref_mobility)
    # Y: the mobility distance
    return Y

def cal_size(Y, gel_conc=1.0, ref_mobility=50, formula='Wu'):
    '''Predict size based on the relative mobility'''
    import math

    gel_para_dict, a, b, c = load_gel_para_dict(gel_conc=gel_conc, formula=formula)

    Y = float(Y)
    gel_conc = float(gel_conc)

    # Y: size (bp)
    # ref_mobility: the mobility distance of the fastest DNA segment

    if gel_conc == 0.5:
	Y = 2.709423 - 0.2691262 * math.log(X+464.4412) 
    elif gel_conc == 1.5:
	Y = 2.322184 - 0.2634794 * math.log(X+48.08739)           
    elif gel_conc == 2.0:
	Y = 2.133323 - 0.2561728 * math.log(X-18.54176)           
    else:
	# Use default gel concentration: 1.0
	Y = 2.397726 - 0.2700602 * math.log(X+73.97887)         
    # X: the mobility distance
    # length: the mobility distance of the fastest DNA segment
    X = X / length
    # here, X was been convert to the relative mobility = mobility distance / length
    Y = math.exp(6.353971 - 1.384176 * math.log(X)) - 474.6539
    # Y: size (bp)
    return Y


def main ():
    '''Test'''
    import sys
    X = 100
    gel_conc = sys.argv[1]
    mobility = cal_mobility(X, gel_conc, ref_mobility=50, formula='Helling')
    print mobility
    mobility = cal_mobility(X, gel_conc, ref_mobility=50)
    print mobility


if __name__ == '__main__':
    main()

