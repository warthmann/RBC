#!/usr/bin/env python
# -*- coding: utf-8 -*- #   
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
                                hid = m.group(4).strip()
                            else:
                                print 'Parsing error!'
                                print line
                                exit()

                        if re.match('^  PPC = ', line):
                            align = line

                            r = re.compile('^  PPC = (.+)%, Size = (\d+) bp, GC content = (.+)%')
                            m = r.search(line)
                            if m:
                                size = m.group(2)
                            else:
                                print 'Parsing error!'
                                print line
                                exit()

                            while 1:
                                line = f.readline()
                                if not line:
                                    break
                                align = align + line
                                if re.search('<<<$', line):
                                    break

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
                        'fid' : fid,
                        'rid' : rid,
                        'hid' : hid,
			'align' : align,
                        'size' : int(size),
                        'desc' : desc,
                        'seq' : seq,
                    }

    return amplicons


def main ():
    filename = parseOptions()
    amplicons = parse(filename)
    for sn in amplicons.keys():
	print amplicons[sn].keys()
        print amplicons[sn]['fid']
        print amplicons[sn]['rid']
        print amplicons[sn]['hid']
        print amplicons[sn]['size']
        print amplicons[sn]['desc']
        print amplicons[sn]['seq']
        print amplicons[sn]['align']

if __name__ == '__main__':
    main()

