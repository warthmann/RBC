#!/usr/bin/env python
# -*- coding: utf-8 -*- #   

import sys
import re
import string
import LightBlastParser as LBP

# Update infomation

# Mar 3, 2009
# Fix the bug which cumulate the Blast record when parsing a set of Blast reports

# May 12, 2008
# Converse Score value (and so on) to number, former is str

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

class nextHsp:
    expect = ''
    bits = ''
    score = ''
    identity = ''
    match = ''
    query_strand = ''
    sbjct_strand = ''
    query_begin = ''
    query_end = ''
    sbjct_begin = ''
    sbjct_end = ''
    qseq = ''
    aseq = ''
    sseq = ''
    sgaps = ''
    qgasp = ''
    mismatch = ''


    def __init__(self, fh, last_line):
        self.__parseHeader(fh, last_line)

    def __parseHeader(self, fh, line):
        hspline = []
        while 1:
            pos = fh.tell()
            if re.match('^\s+Score = ', line):
                r = re.compile('Score =\s+(\S+)\s+bits\s+\((\d+)\)\,\s+Expect\(?\d?\)?\s+=\s+(\S+)')
                m = r.search(line)
                if m:
                    self.bits = float(m.group(1))
                    self.score = int(m.group(2))
                    self.expect = m.group(3)
                    if re.match('^e', self.expect):
                        self.expect = '1' + self.expect
                        self.expect = float(self.expect)
                else:
                    print 'Parsing error in %s' % line
                    exit()

            if re.match('^\s+Identities = ', line):
                r = re.compile('Identities =\s+(\d+)/(\d+) ')
                m = r.search(line)
                if m:
                    self.identity = int(m.group(1))
                    self.match = int(m.group(2))
                else:
                    print 'Parsing error in %s' % line
                    exit()

            if re.match('^\s+Strand = ', line):
                r = re.compile('Strand =\s+(\S+)\s+\/\s+(\S+)')
                m = r.search(line)
                if m:
                    self.query_strand = m.group(1)
                    self.sbjct_strand = m.group(2)
                else:
                    print 'Parsing error in %s' % line
                    exit()

            if re.match('^Query: ', line):
                hspline.append(line) # Query line
                ali_line = fh.readline()  # Would be the alignment line
                if re.match('^Sbjct: ', ali_line): 
                    hspline.append("") # dummy line, this is a -noseq option
                    hspline.append(ali_line) # so store a fake alignment and real sbjct
                else:
                    hspline.append(ali_line)
                    sbj_line = fh.readline() # Sbjct line
                    hspline.append(sbj_line)


            line = fh.readline()
            if not line:
                break
            if re.match('^\s+Score =\s+', line) and hspline:
                self.parseAlignment(hspline)
                fh.seek(pos)
                break
            if re.match('^>\S+', line) and hspline:
                self.parseAlignment(hspline)
                fh.seek(pos)
                break
            if re.match('^BLAST?', line) and hspline:
                self.parseAlignment(hspline)
                fh.seek(pos)
                break
            if re.match('^  Database:', line) and hspline:
                self.parseAlignment(hspline)
                break

    def parseAlignment(self, hspline):
        qb = -1
        qe = -1
        sb = -1
        se = -1
        query_seq = []
        align_seq = []
        sbjct_seq = []
        th = 0
        while th < len(hspline):
            line = hspline[th]
            r = re.compile('^Query:\s+(\d+)\s+(\S+)\s+(\d+)')
            m = r.search(line)
            if not m:
                print "Alignment lines parsing error!\n %s" % line
                exit()
            if qb < 0:
                qb = int(m.group(1))
            qseq = m.group(2)
            qe = int(m.group(3))

            offset = line.index(qseq)
            if hspline[th + 1]: # Alignment line
                aseq = hspline[th + 1][offset : (offset + len(qseq))]
            else:
                print 'Parsing Error. \n %s' % hspline[th+1]
                exit()
            line = hspline[th + 2]
            r = re.compile('^Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)')  # Sbjct line
            m = r.search(line)
            if not m:
                print "Alignment lines parsing error!\n %s" % line
                exit()
            if sb < 0:
                sb = int(m.group(1))
            sseq = m.group(2)
            se = int(m.group(3))

            query_seq.append(qseq)
            align_seq.append(aseq)
            sbjct_seq.append(sseq)

            th = th + 3

        # the HSP variables
        if qb < 0 or qe < 0 or se < 0 or sb < 0:
            raise ParsingError, "Parsing Error in Hsp alignment"
        self.query_begin = qb
        self.query_end = qe
        self.sbjct_begin = sb
        self.sbjct_end = se
        self.qseq = ''.join(query_seq).upper()
        self.aseq = ''.join(align_seq)
        self.sseq = ''.join(sbjct_seq).upper()
        self.qgaps = qseq.count('-')
        self.sgaps = sseq.count('-')
        self.mismatch = aseq.count(' ')


class nextHit:
    hsps = []
    hit_id = ''
    hit_desc = ''
    hit_length = ''
    def __init__(self, fh, last_line):
        self.hsps = self.__parseHeader(fh, last_line)

    def __parseHeader(self, fh, line):
        tmp_hsps = []
        while 1:
            if re.match('^>\S+', line):
                hit_full = re.sub('>', '', line)
                self.hit_id = hit_full.split(' ')[0]
                if len(line.split()) > 1:
                    hit_full = ' '.join(hit_full.split(' ')[1:])
                else:
                    hit_full = ''
                while 1:
                    line = fh.readline()
                    r = re.compile('^\s+Length =\s+(\d+)')
                    m = r.search(line)
                    if m:
                        self.hit_length = m.group(1)
                        break # Smart
                    else:
                        line = line.strip()
                        hit_full = hit_full + line
                self.hit_desc = ' '.join(hit_full.split('\n')[:])

            if re.match('^\s+Score =\s+', line):
                last_line = line
                tmp_hsps.append(nextHsp(fh, last_line))

            pos = fh.tell()
            line = fh.readline()
            if not line:
                break
            if re.match('^>\S+', line):
                fh.seek(pos)
                break
            if re.match('^BLAST?', line):
                fh.seek(pos)
                break
            if re.match('^  Database:', line):
                break
        return tmp_hsps


class nextRecord:
    hits = []
    version = ''
    reference = ''
    database = ''
    db_sequences = ''
    db_letters = ''
    query_id = ''
    query_desc = ''
    query_letters = ''

    def __init__(self, fh, last_line):
        self.hits = self.__parseHeader(fh, last_line)

    def __parseHeader(self, fh, line):
        tmp_hits = []
        while 1:
            if re.match('^BLAST?', line):
                    self.version = line.rstrip() # version
            # Parsing for reference
            if re.match("^Reference: ", line):
                reference = line.split(' ')[1:]  # reference
                reference = ' '.join(reference)
                while 1:
                    line = fh.readline()
                    if re.match('\n', line):
                        # After the reference, there must be a space line exists
                        break
                    else:
                        reference = reference + line
                self.reference = reference
            # Parsing for query information
            if re.match("^Query= ", line):
                query_full = line
                while 1:
                    line = fh.readline()
                    r = re.compile('^\s+\((\S+)\s+letters\)$')
                    m = r.search(line)
                    if m:
                        self.query_letters = m.group(1)
                        self.query_letters = re.sub('[^\d]+', '', self.query_letters)
                        self.query_letters = int(self.query_letters)
                        break # Smart
                    else:
                        query_full = query_full + line
                self.query_id = query_full.split()[1]
                if len(query_full.split()) > 1:
                    query_desc = ' '.join(query_full.split()[2:])
                else:
                    query_desc = ''
                self.query_desc = query_desc

            # 
            # Parsing for database infomation
            if re.match("^Database: ", line):
                database_full = line
                while 1:
                    line = fh.readline()
                    r = re.compile('^\s+(.+) sequences; (.+) total letters')
                    m = r.search(line)
                    if m:
                        self.db_sequences = m.group(1)
                        self.db_letters = m.group(2)
                        break # Smart
                    else:
                        database_full = database_full + line
                self.database = ' '.join(database_full.split()[1:])

            # Parsing for 'No hits found'
            if re.search("No hits found", line):
                self.hits = []
            
            if re.match("^>\S+", line):
                last_line = line
                tmp_hits.append(nextHit(fh, last_line))

            pos = fh.tell()
            line = fh.readline()
            if not line:
                break
            if re.match('^BLAST?', line):
                fh.seek(pos)
                break
            if re.match('^  Database:', line):
                break
        return tmp_hits



class results:
    records = []
    def __init__(self, fh):
        self.records = []
        self.__nextRecord(fh)

    def __nextRecord(self, fh):
        while 1:
            line = fh.readline()
            if not line:
                break
            if re.match('^BLAST?', line):
                last_line = line
                self.records.append(nextRecord(fh, last_line))
            if re.match('^  Database:', line):
                break

def parse(fh):
    return results(fh)

def main ():
    # Usage
    blast_file = parseOptions()
    fh = open(blast_file, 'r')
    result = LBP.parse(fh)
    fh.close()

    for record in result.records:
        print record.version
        print record.reference
        print record.query_id, record.query_desc
        print '\t', record.query_letters
        print '\n'
        print record.database
        print record.db_sequences, record.db_letters

        for hit in record.hits:
            print '>%s %s' % (hit.hit_id, hit.hit_desc)
            print '\tLength = %s' % hit.hit_length
            print '\n'
            for hsp in hit.hsps:
                print ' Score = %s bits (%s), Expect = %s' % (hsp.bits, hsp.score, hsp.expect)
                print ' Identities = %s/%s' % (hsp.identity, hsp.match)
                print ' Strand = %s / %s'  % (hsp.query_strand, hsp.sbjct_strand)
                print ' qb = %s, qe = %s; sb = %s, se = %s' % (hsp.query_begin, hsp.query_end, hsp.sbjct_begin, \
                                                               hsp.sbjct_end)
                print '\n'
                print hsp.qseq
                print hsp.aseq
                print hsp.sseq
                print '\n'
                print '\n'


if __name__ == '__main__':
    main()
