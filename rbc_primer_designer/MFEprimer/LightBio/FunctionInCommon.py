#!/usr/bin/env python
# -*- coding:utf-8 -*-

#-----------------------------------------------------
# File: function.py
# Date: 2008-09-01
# Description: Collect some useful function here
#-----------------------------------------------------

__version__ = '1.0'
__author__  = 'Wubin Qu <quwubin@gmail.com> @ZCGLAB @BMI @CHINA'

def print_seq(seq, width=80):
    '''Print biological sequence with constant width, such as 80'''
    import os

    line = ''
    start = 0
    while start <= len(seq):
        seq_line = seq[start : (start + width)]
        start = start + width
        line = line + seq_line + os.linesep

    return line 

def num_group(number):
    '''Make 123456 like 123,456'''
    s = '%d' % number
    groups = []
    while s and s[-1].isdigit():
        groups.append(s[-3:])
        s = s[:-3]
        
    return s + ','.join(reversed(groups))

def random_string(length):
    import string
    import random
    rstr = ''.join(random.choice(string.letters) for i in xrange(length))
    return rstr

def get_size_range(size):
    Y = cal_mobility(size)
    # Set 2 mm as the distance which the bands can be \
            #seperated by naked eyes
    Xmin = cal_size(Y+2)
    Xmax = cal_size(Y-2)

    return Xmin, Xmax

def cal_mobility(X, length=50):
    import math

    X = float(X)
    # X: size (bp)
    # length: the mobility distance of the fastest DNA segment
    Y = math.exp(4.606033 - 0.7245908 * math.log(X + 474.6539))
    # Y: the relative mobility = mobility distance / length
    Y = Y * length
    # Y: the mobility distance
    return Y

def cal_size(X, length=50):
    import math

    X = float(X)
    # X: the mobility distance
    # length: the mobility distance of the fastest DNA segment
    X = X / length
    # here, X was been convert to the relative mobility = mobility distance / length
    Y = math.exp(6.353971 - 1.384176 * math.log(X)) - 474.6539
    # Y: size (bp)
    return Y

def codec_read_file(file_name, codec_char = 'utf-8'):
    import codecs

    fh = codecs.open(file_name, encoding=codec_char)
    file_content = fh.read()
    fh.close()
    return file_content


def read_from_file(file_name):
    fh = open(file_name)
    file_content = fh.read()
    fh.close()
    return file_content

def write_to_file(file_name, file_content, model):
    fh = open(file_name, model)
    fh.write(file_content)
    fh.close()

def connect_mysql(host, user, passwd, dbname):
    ''' Connect to the MySQL database '''
    import MySQLdb
    db = MySQLdb.Connection(host, user, passwd, dbname)
    cur = db.cursor()

    return cur, db

def debug(s):
    fh = open('debug.tmp', 'a')
    fh.write(str(s))
    fh.write('\n')
    fh.close()

def close_mysql(cur, db):
    ''' Close the database '''
    cur.close()
    db.close()

def parse_csv(filename, delimiter_char=','):
    import csv

    csvfile = open(filename)
    reader = csv.reader(csvfile, delimiter=delimiter_char)
    has_header = True
    dicts = []
    for fields in reader:
        if has_header:
            header = fields
            has_header = False
        else:
            dicts.append({})
            for i, f in enumerate(fields):
                dicts[-1][header[i]] = f
    return dicts

def set_cache(content, path, file_name):
    import os
    import shelve

    full_path = path + os.sep + file_name
    d = shelve.open(full_path)
    d[file_name] = content
    d.close()

def get_cache(path, file_name):
    import os
    import shelve

    full_path = os.path.join(path, file_name)
    d = shelve.open(full_path)
    content = d[file_name]
    d.close()

    return content

def main():
    '''Test'''
    min_size, max_size = get_size_range(100)
    print min_size, max_size

if __name__ == '__main__':
    main()
