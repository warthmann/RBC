#!/usr/bin/env python
# -*- coding: utf-8 -*- #   

# Parsing FASTA alignment (-m 2) output for ThermoFASTA
# Note: this parser may not work for other applications

# Author: Wubin Qu [quwubin@gmail.com]
# Date: 2009.09.16 16:54:58
# License: GPL(v3)

import sys
import re

class nextHsp:
    hit_id = ''
    hit_desc = ''
    qb = ''
    qe = ''
    sb = ''
    se = ''
    qseq = ''
    sseq = ''

    def __init__(self, fh, last_line, query_id):
        self.__parseHeader(fh, last_line, query_id)

    def __parseHeader(self, fh, line, query_id):
        while 1:
            pos = fh.tell()
            if re.match('^>>\S+', line):
                tmp_row = line.split()
                self.hit_id = tmp_row[0]
                self.hit_id = re.sub('>>', '', self.hit_id)
                try:
                    self.hit_desc = ' '.join(tmp_row[1:-3])
                except:
                    pass

            if re.match('^banded Smith-Waterman', line):
                pos_str = line.split()[-1]
                pos_str = re.sub('\(|\)', '', pos_str)
                self.qb, self.qe = pos_str.split(':')[0].split('-')
                self.sb, self.se = pos_str.split(':')[1].split('-')
                self.qb = int(self.qb)
                self.qe = int(self.qe)
                self.sb = int(self.sb)
                self.se = int(self.se)
            
            if line.startswith(query_id):
                self.qseq = line.split()[1]

            try:
                if line.startswith(self.hit_id[:19]):
                    self.sseq = line.split()[1]
            except:
                print 'sseq parsing error:\n%s' % line
                exit()

            line = fh.readline()
            if not line:
                break

            if re.match('^>>\S+', line):
                fh.seek(pos)
                break
            if re.match('^\d+\s+residues\s+in\s+\d+\s+library\s+sequences', line):
                fh.seek(pos)
                break
            if re.match('^\s+\d+>>>', line):
                fh.seek(pos)
                break

        assert self.hit_id and self.qb and self.qe and self.sb and self.se and self.qseq and self.sseq, 'FASTA report query line parsing error'

class nextRecord:
    '''Parsing (Spliting) the file based on the different query'''
    hsps = []
    query_sn = ''
    query_id = ''
    query_desc = ''
    query_letters = ''

    def __init__(self, fh, last_line):
        self.hsps = self.__parseHeader(fh, last_line)

    def __parseHeader(self, fh, line):
        tmp_hsps = []
        while 1:
            if re.match('^\s+\d+>>>', line):
                r = re.compile('^\s+(\d+)>>>(.+)-\s+(\d+)\s+nt')
                m = r.search(line)
                try:
                    self.query_sn = int(m.group(1))
                    self.query_id = m.group(2).split()[0].strip()
                    self.query_letters = int(m.group(3))
                    try:
                        self.query_desc = ' '.join(m.group(2).split()[1:])
                    except:
                        pass
                except:
                    print 'Parsing error in line: \n%s\n' % line
                    exit()

            if re.match('^>>\S+', line):
                last_line = line
                tmp_hsps.append(nextHsp(fh, last_line, self.query_id))

            pos = fh.tell()
            line = fh.readline()
            if not line:
                break
            if re.match('^\d+\s+residues\s+in\s+\d+\s+library\s+sequences', line):
                fh.seek(pos)
                break
            if re.match('^\s+\d+>>>', line):
                fh.seek(pos)
                break

        assert self.query_sn and self.query_id and self.query_letters, 'FASTA report query line parsing error'

        return tmp_hsps

class results:
    '''Begin to parse'''
    records = []
    version = ''
    reference = 'W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448'
    db_letters = ''
    db_sequences = ''

    def __init__(self, fh):
        self.records = []
        self.__nextRecord(fh)

    def __nextRecord(self, fh):
        while 1:
            line = fh.readline()
            if not line:
                break
            if re.match('\s+version', line):
                self.version = ' '.join(line.split()[1:])

            if re.match('^\s+\d+>>>', line):
                last_line = line
                self.records.append(nextRecord(fh, last_line))

            if re.match('^\d+\s+residues\s+in\s+\d+\s+library\s+sequences', line):
                # The end of the file: summary inforation
                tmp_row = line.split()
                self.db_letters = int(tmp_row[0])
                self.db_sequences = int(tmp_row[3])
        
        assert self.version and self.db_letters and self.db_sequences, 'FASTA report format error'
    
def parse(fh):
    '''Parse the file and return the results'''
    return results(fh)

def main ():
    '''Test the parser'''
    results = parse(open(sys.argv[1]))
    print results.version
    print results.reference
    print results.db_letters
    print results.db_sequences
    for record in results.records:
        print '\tquery_sn', record.query_sn
        print '\tquery_id', record.query_id
        print '\tquery_desc', record.query_desc
        print '\tquery_letters', record.query_letters
        for hsp in record.hsps:
            print '\t\thit_id:', hsp.hit_id
            print '\t\thit_desc:', hsp.hit_desc
            print '\t\tqb:', hsp.qb
            print '\t\tqe:', hsp.qe
            print '\t\tsb:', hsp.sb
            print '\t\tse:', hsp.se
            print '\t\tqseq:', hsp.qseq
            print '\t\tsseq:', hsp.sseq

if __name__ == '__main__':
    main()

