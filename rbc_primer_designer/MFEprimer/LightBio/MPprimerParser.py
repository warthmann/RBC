#!/usr/bin/env python
# -*- coding: utf-8 -*- #   

Author = 'Wubin Qu <quwubin@gmail.com>, BIRM, China'
Date = '2009-11-23'
License = 'GPL v3'
Version = '1.0'

import sys
import re
import string

def parseOptions ():
    if len(sys.argv) < 2:
	print """
	USAGE:
	python program.py file
	"""
	sys.exit(1)
    else:
	filename = sys.argv[1]
    return filename

def parse(filename):
    psc = {}
    
    f = open(filename, 'r')
    while 1:
        line = f.readline()
        if not line:
            break
        line = line.rstrip()
        if re.match('^Details of', line):
            while 1:
                line = f.readline()
                if not line:
                    break
                line = line.rstrip()
                next_psc = 1
                t_dict = {}
                psc_sn = None
                if next_psc:
                    next_psc = 0
                    end_psc = 0
                    while not end_psc: 
                        line = f.readline()
                        if not line:
                            break
                        line = line.rstrip()
                        if re.match('^PSC\s+\d+:', line):
                            r = re.compile('^PSC\s+(\d+):')
                            m = r.search(line)
                            if m:
                                psc_sn = m.group(1)
                            else:
                                print 'Parsing error!'
                                print line
                                exit()

                        if re.match('^T\d+', line):
                            while 1:
                                r = re.compile('^T(\d+):\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+/\s+(\S+)\s+(\S+)\s+/\s+(\S+)')
                                m = r.search(line)
                                if m:
                                    t_sn = m.group(1)
                                    t_id = m.group(2)
                                    penalty = m.group(3)
                                    size = m.group(4)
                                    fp_tm = m.group(5)
                                    rp_tm = m.group(6)
                                    fp = m.group(7)
                                    rp = m.group(8)
                                else:
                                    print 'Parsing error!'
                                    print line
                                    exit()

                                t_sn = int(t_sn)
                                t_dict[t_sn] = {
                                    't_id' : t_id,
                                    'penalty' : penalty,
                                    'size' : int(size),
                                    'fp_tm' : float(fp_tm),
                                    'rp_tm' : float(rp_tm),
                                    'fp' : fp,
                                    'rp' : rp,
                                }

                                line = f.readline()
                                line = line.rstrip()
                                if not line:
                                    break

                        if re.match('^\*+', line):
                            end_psc = 1
                            next_psc = 1

                    if next_psc and psc_sn:
                        psc[int(psc_sn)] = t_dict
                        psc_sn = None

    return psc


def main ():
    filename = parseOptions()
    psc = parse(filename)
    for psc_sn in psc.keys():
        print 'PSC: %s' % psc_sn
        for t_sn in psc[psc_sn].keys():
            print '\tT: ', t_sn
            print '\t', psc[psc_sn][t_sn]['t_id']
            print '\t', psc[psc_sn][t_sn]['penalty']
            print '\t', psc[psc_sn][t_sn]['size']
            print '\t', psc[psc_sn][t_sn]['fp_tm']
            print '\t', psc[psc_sn][t_sn]['rp_tm']
            print '\t', psc[psc_sn][t_sn]['fp']
            print '\t', psc[psc_sn][t_sn]['rp']

if __name__ == '__main__':
    main()

