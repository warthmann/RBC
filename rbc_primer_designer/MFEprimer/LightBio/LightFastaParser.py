#!/usr/bin/env python
__version__ = '1.0'
__date__ = 'May 18, 2008'
__author__ = 'Wubin Qu <quwubin@gmail.com> @ZCGLAB @BMI'
__url__ = 'http://quwubin.wordpress.com.cn/'
__license__ = 'GPL'

import sys
import re
import LightFastaParser as LFP

class FastaRecord(object):
    def __init__(self, fh, last_line):
        self.__parse(fh, last_line)

    def __parse(self, fh, line):
        sequences = []
        while 1:
            if re.match('^>\S+', line):
                self.title = line[1:]
                if len(self.title.split()) > 1:
                    self.id = self.title.split()[0]
                    self.desc = self.title.split()[1:]
                else:
                    self.id = self.title
                    self.desc = ''

            pos = fh.tell()
            line = fh.readline()
            line = line.strip()
            if not line:
                break
            if re.match('^>\S+', line):
                fh.seek(pos)
                break
            sequences.append(line)
        self.sequence = ''.join(sequences)
        self.length = len(self.sequence)

def parse(fh):
    records = []
    while 1:
        line = fh.readline()
        line = line.strip()
        if not line:
            break
        if re.match('^>\S+', line):
            last_line = line
            records.append(FastaRecord(fh, last_line))
    return records


# The master test function
def test():
    filename = sys.argv[1]
    fh = open(filename, 'r')
    records = LFP.parse(fh)
    fh.close()
    for record in records:
        #print record.id
        print '>%s' % record.title
        print record.sequence
        #print record.length
    
if __name__ == "__main__":
    test()
