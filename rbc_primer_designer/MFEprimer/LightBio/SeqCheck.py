#!/usr/bin/python
'''A module for validating the format of biological sequences, such 
as Fasta format, etc.'''

Author = 'Wubin Qu <quwubin@gmail.com'
Date = '2009-10-19'
License = 'GPL v3'

import sys, re
import math

def print2err(msg):
    import sys
    print >> sys.stderr, msg

def FastaCheck(seq):
    '''Validating the Fasta format of given sequences'''
    import sys, re, math

    seq = seq.strip()

    if not re.match('>', seq):
        print2err("Please input the sequence in FASTA format")
        exit()

    if not re.search('[atcgnATCGN]+', seq):
        print2err("Please input nucleotide sequence in FASTA format")
        exit()

    lines = seq.split('\n')
    for line in lines:
	if line.strip() == '':
	    lines.remove(line)
    seq = '\n'.join(lines)

    fastas = seq.split('>')
    line_list = []
    for fasta in fastas:
        if fasta.strip():
            lines = fasta.split('\n')
            if lines[0].strip():
                desc_line = '>' + lines[0].strip()
                line_list.append(desc_line)
            seq_lines = ''
            if len(lines[1:]) > 0:
                for line in lines[1:]:
                    line = line.strip()
                    seq_lines = seq_lines + line

            if seq_lines:
                line_list.append(seq_lines)

    if math.fmod(len(line_list), 2):
        print2err("Please input the sequence in FASTA format")
        exit()
    else:
        count = 0
        for line in line_list:
            count = count + 1
            if math.fmod(count, 2): # > line 
                if not re.match('>', line):
                    print2err("Please input the sequence in FASTA format")
                    exit()
            else:  # sequence line
                if re.search('[^atcgnATCGN]+', line):
                    print2err("Please input nucleotide sequence in FASTA format")
                    exit()

    return seq

def main():
    '''Test'''
    fh = open(sys.argv[1])
    seq = fh.read()
    seq = FastaCheck(seq)
    print seq
    fh.close()

if __name__ == '__main__':
    main()
