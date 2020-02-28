#!/usr/bin/env python
# -*- coding: utf-8 -*- #   
'''A parser for FASTA m 9 output report'''

Author = 'Wubin Qu <quwubin@gmail.com>, BIRM, China'
Date = '2009-12-10'
License = 'GPL v3'
Version = '1.0'

import sys, re, os

def parse(out_text):
    '''Parse FASTA report in memory
    Return a Dict 
    '''
    parts = re.split('\d+>>>', out_text)
    rs = {}  # r for fasta reports
    rs['reference'] = 'W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448'

    for line in parts[0].split('\n'):
        if re.search('^\s+version', line):
            version = ' '.join(line.strip().split()[1:])
            rs['version'] = version


    rs['r'] = {}
    for part in parts[1:]:
        flag = False
        hit_line = []
        lines = part.split('\n')
        qid = lines[0].split()[0] # Risk
        rs['r'][qid] = {}
        for line in lines:
            if re.search('^The best scores are', line):
                flag = True
                continue

            if flag:
                if not line.strip():
                    break
                hit_line.append(line)
        if flag:
            hits = []
            for line in hit_line:
                items = line.split()
                hid = items[0]
                qb = items[-11]
                qe = items[-10]
                sb = items[-7]
                se = items[-6]
                hit_length = items[-4]
                orientation = items[-19][1]
                score = items[-17]
                evalue = items[-16]

                tmp = {}
                tmp = {
                    'hid' : hid,
                    'qb' : int(qb),
                    'qe' : int(qe),
                    'sb' : int(sb),
                    'se' : int(se),
                    'ot' : orientation,
                    'score' : score,
                    'evalue' : evalue,
                    'hit_length' : int(hit_length),
                }
                hits.append(tmp)

            rs['r'][qid] = hits

    return rs

def main ():
    '''Main test function'''
    rs = parse(open(sys.argv[1]).read())
    print rs['version']
    print rs['reference']

if __name__ == '__main__':
    main()

