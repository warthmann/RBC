#!/usr/bin/env python
# -*- coding: utf-8 -*- #   
'''Generate permutations with repeat items in an array'''

Author = 'Wubin Qu <quwubin@gmail.com>, BIRM, China'
Date = '2010-1-13'
License = 'GPL v3'
Version = '1.0'

import sys, re, os
array_out = []

def write(rec, N, array):
    '''Output'''
    out = []
    for i in range(N):
        out.append(array[rec[i]])

    array_out.append(out)

def arrange(rec, used, depth, N, array):
    '''Arrange'''
    if depth >= N:
        write(rec, N, array)
    else:
        found_num = sys.maxint
        for i in range(N):
            if used[i] == 0 and array[i] < found_num:
                rec[depth] = i
                found_num = array[i]
                used[i] = 1
                arrange(rec, used, depth+1, N, array)
                used[i] = 0

def main ():
    '''Main'''
    test_array = [1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    #test_array = [1, 1, 2, 2]
    #test_array = [2, 2, 1, 1]
    #test_array = [1, 1, 2]
    #test_array = [1, 0]
    #test_array = [1, 0]
    #test_array = [2, 2, 1, 1, 1]
    #test_array = [1, 1, 0, 0]
    test_array = [1, 1, 1, 1]
    N = len(test_array)
    count = 0
    rec = [0 for col in range(N+1)]
    used = [0 for col in range(N+1)]
    depth = 0
    arrange(rec, used, depth, N, test_array)
    print len(array_out)
    for i in array_out:
        print i


if __name__ == '__main__':
    main()

