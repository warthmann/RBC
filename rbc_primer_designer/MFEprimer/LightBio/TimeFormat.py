#!/usr/bin/env python
# -*- coding: utf-8 -*- #   
'''A module to manipulate time and format output'''

# Author: Wubin Qu [quwubin@gmail.com]
# Date: 2009-10-09

def seconds_format(s):
    '''Format output senconds as Hour:Minutes:Seconds'''
    if s/60 != 0:
        min = s/60
        sec = s%60
    else:
        min = 0
        sec = s
    return (min, sec)

def main ():
    import sys
    m, s = seconds_format(int(sys.argv[1]))
    print '%s m %s s' % (m, s)


if __name__ == '__main__':
    main()

