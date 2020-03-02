#!/usr/bin/env python
# -*- coding: utf-8 -*- #   
'''A parser for Primer3 output

Support File and String (passed to parse function) parameters.

Author: Wubin Qu <quwubin@gmail.com>
Date: 2009-09-28
'''

import sys
import re
from types import *

def parse(output):
    '''Parse function for Primer3 output'''
    p3out = []
    result_string = ''
    if type(output) is FileType:
        result_string = output.read()
    elif type(output) is StringType:
        result_string = output
    else:
        print 'Not legal Primer3 output found'
        exit()

    if not result_string:
        print 'Not legal Primer3 output found'
        exit()

    for report in result_string.split('\n=\n'):

        if not report:
            continue

        item_dict = {}
        for line in report.split('\n'):
            if not line.strip():
                continue
            key, value = line.split('=')
            item_dict[key.strip()] = value.strip()
        p3out.append(item_dict) 

    return p3out

def main ():
    '''Test function'''
    p3out = parse(open(sys.argv[1]))
    for report in p3out:
        ok = report['PRIMER_PAIR_EXPLAIN']
        print ok
        for key, value in report.items():
            if key == 'PRIMER_PAIR_EXPLAIN':
                print '***************************'

            if re.search('EXPLAIN', key):
                pass
                #print key, value


if __name__ == '__main__':
    main()

