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
                        'desc' : desc,
                        'seq' : seq,
                    }

    return amplicons


def main ():
    filename = parseOptions()
    amplicons = MP.parse(filename)
    for sn in amplicons.keys():
        print amplicons[sn]['desc']
        print amplicons[sn]['seq']

if __name__ == '__main__':
    main()

