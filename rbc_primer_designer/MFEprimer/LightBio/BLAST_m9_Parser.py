#!/usr/bin/env python
# -*- coding: utf-8 -*- #   
'''A parser for BLAST tabular with comment lines output'''

Author = 'Wubin Qu <quwubin@gmail.com>, BIRM, China'
Date = '2009-12-18'
License = 'GPL v3'
Version = '1.0'

import sys, re, os

def parse(output):
    '''Parse the tabular with comment lines (m 9) output'''
    rs = {}
    rs['reference'] = '''Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs",  Nucleic Acids Res. 25:3389-3402.'''
    rs['r'] = {}
    for line in output.split('\n'):
        if not rs.has_key('version'):
            if re.match('^#\s+BLAST', line):
                rs['version'] = ' '.join(line.split()[2:])
        if re.match('^#', line):
            continue
        if not line.strip():
            continue
        items = line.split()
        qid = items[0]
        hid = items[1]
        qb = int(items[6])
        qe = int(items[7])
        sb = int(items[8])
        se = int(items[9])
        score = items[-1]
        evalue = items[-2]

        if se > sb:
            ot = 'f'
        else:
            ot = 'r'

        tmp = {
            'hid' : hid,
            'qb' : qb,
            'qe' : qe,
            'sb' : sb,
            'se' : se,
            'score' : score,
            'evalue' : evalue,
            'hit_length' : None,
            'ot' : ot,
        }

        if not rs['r'].has_key(qid):
            rs['r'][qid] = []
            rs['r'][qid].append(tmp)
        else:
            rs['r'][qid].append(tmp)

    return rs

def main ():
    '''Main'''
    rs = parse(open(sys.argv[1]).read())
    print rs['version']
    print rs['reference']
    print rs['r'].keys()
    for qid in rs['r'].keys():
        print rs['r'][qid]['evalue']
        print rs['r'][qid]['ot']
        print rs['r'][qid]['sb']
        print rs['r'][qid]['se']



if __name__ == '__main__':
    main()

