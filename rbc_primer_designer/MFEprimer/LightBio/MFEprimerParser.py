#!/usr/bin/env python
# -*- coding: utf-8 -*- #   
import sys
import re
import string
import MFEprimerParser as MP

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
    amplicons = {}
    
    f = open(filename, 'r')
    while 1:
        line = f.readline()
        if not line:
            break
        line = line.rstrip()
        if re.match('^Details for', line):
            while 1:
                line = f.readline()
                if not line:
                    break
                line = line.rstrip()
                next_amplicon = 1
                if next_amplicon:
                    next_amplicon = 0
                    end_amplicon = 0
                    while not end_amplicon: 
                        line = f.readline()
                        if not line:
                            break
                        line = line.rstrip()
                        if re.match('^\d+: ', line):
                            r = re.compile('^(\d+): (.+) \+ (.+) ==> (.+)')
                            m = r.search(line)
                            if m:
                                sn = m.group(1)
                                fid = m.group(2)
                                rid = m.group(3)
                                hid = m.group(4)
                            else:
                                print 'Parsing error!'
                                print line
                                exit()

                        if re.match('^  PPC = ', line):
                            r = re.compile('^  PPC = (.+)%, Size = (\d+) bp, GC content = (.+)%')
                            m = r.search(line)
                            if m:
                                ppc = m.group(1)
                                size = m.group(2)
                                gc_content = m.group(3)
                            else:
                                print 'Parsing error!'
                                print line
                                exit()

                        if re.match('^  FP: ', line):
                            r = re.compile('^  FP: 3\'ΔG = (.+)\(kcal/mol\), Tm = (.+) \(°C\)')
                            m = r.search(line)
                            if m:
                                fp_deltaG = m.group(1)
                                fp_tm = m.group(2)
                                if fp_tm == '-':
                                    fp_tm = 0
                            else:
                                print 'Parsing error!'
                                print line
                                exit()

                        if re.match('^  RP: ', line):
                            r = re.compile('^  RP: 3\'ΔG = (.+)\(kcal/mol\), Tm = (.+) \(°C\)')
                            m = r.search(line)
                            if m:
                                rp_deltaG = m.group(1)
                                rp_tm = m.group(2)
                                if rp_tm == '-':
                                    rp_tm = 0
                            else:
                                print 'Parsing error!'
                                print line
                                exit()

                        if re.match('^  Binding sites: ', line):
                            r = re.compile('^  Binding sites: (\d+)\(.+\) ... (\d+)\(.+\)')
                            m = r.search(line)
                            if m:
                                start = m.group(1)
                                stop = m.group(2)
                            else:
                                print 'Parsing error!'
                                print line
                                exit()

                        if re.match('^>\d.+', line):
                            desc = line
                            seq = ''
                            while 1:
                                line = f.readline()
                                if not line:
                                    break
                                if line == '\n':
                                    end_amplicon = 1
                                    next_amplicon = 1
                                    break
                                seq = seq + line

                if next_amplicon:
                    tmp = int(sn) - 1
                    amplicons[tmp] = {
                        'fp_deltaG' : float(fp_deltaG),
                        'rp_deltaG' : float(rp_deltaG),
                        'fp_tm' : float(fp_tm),
                        'rp_tm' : float(rp_tm),
                        'fid' : fid,
                        'start' : int(start),
                        'stop' : int(stop),
                        'forward' : fid,
                        'rid' : rid,
                        'reverse' : rid,
                        'hid' : hid,
                        'ppc' : ppc,
                        'ppmi' : ppc,
                        'size' : int(size),
                        'length' : int(size),
                        'gc_content' : float(gc_content),
                        'desc' : desc,
                        'seq' : seq,
                    }

    return amplicons


def main ():
    filename = parseOptions()
    amplicons = MP.parse(filename)
    for sn in amplicons.keys():
        print amplicons[sn]['ppc']
        print amplicons[sn]['fp_deltaG']
        print amplicons[sn]['rp_deltaG']
        print amplicons[sn]['fp_tm']
        print amplicons[sn]['rp_tm']
        print amplicons[sn]['fid']
        print amplicons[sn]['start']
        print amplicons[sn]['stop']
        print amplicons[sn]['forward']
        print amplicons[sn]['rid']
        print amplicons[sn]['reverse']
        print amplicons[sn]['hid']
        print amplicons[sn]['ppmi']
        print amplicons[sn]['size']
        print amplicons[sn]['length']
        print amplicons[sn]['gc_content']
        print amplicons[sn]['desc']
        print amplicons[sn]['seq']

if __name__ == '__main__':
    main()

