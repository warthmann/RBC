#!/usr/bin/env python
# -*- coding:utf-8 -*-

#-----------------------------------------------------
# File: ListParser.py
# Date: 2008-06-30
# Description: 
#-----------------------------------------------------

__version__ = '1.0'
__author__  = 'Wubin Qu <quwubin@gmail.com> @ZCGLAB @BMI @CHINA'
__blog__    = 'http://quwubin.cnblogs.com'
__license__ = 'GPL'

import sys
import re
import ListParser as LP

class info(object):
    def __init__(self, headers, contents):
        self.__parse(headers, contents)

    def __parse(self, headers, contents):
        items = {}
        for i in range(len(headers)):
            items[headers[i]] = contents[i]
        self.items = items

def parse(fh):
    contents = []
    headers = []
    while 1:
        line = fh.readline()
        line = line.strip()
        if not line:
            break

        if re.match('^#', line):
            line = re.sub('^#', '', line)
            headers = line.split('\t')

            for i in range(len(headers)):
                contents.append([])
        else:
            rows = line.split('\t')
            for i in range(len(headers)):
                contents[i].append(rows[i])

    return info(headers, contents)

# The master test function
def test():
    filename = sys.argv[1]
    fh = open(filename, 'r')
    result = LP.parse(fh)
    fh.close()
    print result.items

if __name__ == "__main__":
    test()
    print "Done."
