#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

"""
    MFEprimer.py
    Dec 4, 2008
    ~~~~~~~~

    The MFEprimer.py program is a standalone version of the MFEprimer web 
    interface on http://biocompute.bmi.ac.cn/MFEprimer/.
    Copyright @ 2009. Wubin Qu & Chenggang Zhang. All Rights Reserved.
"""


__version__ = '1.7'
__date__ = 'Jan-28-2010'
__author__ = 'Wubin Qu <quwubin@gmail.com> @CZlab @BMI @CHINA'
__license__ = 'GPL'
__reference__ = '''Wubin Qu, Zhiyong Shen, Dongsheng Zhao, Yi Yang, Chenggang Zhang. (2009)
MFEprimer: Multiple factor evaluation of the specificity of PCR primers, 
Bioinformatics 25(2), 276-278.
'''

import sys
import subprocess
import os
import getopt
import string
import re
import math
src_path = os.path.split(os.path.realpath(sys.argv[0]))[0]
blast_home = src_path + '/bin'
sys.path.append(src_path)
import LightBio.LightFastaParser as LFP
import LightBio.LightBlastParser as LBP
import LightBio.LightSeq as LS
import LightBio.ThermodynamicsParameters as TP
import LightBio.TmDeltaG as TD
import time
import signal
import textwrap
import LightBio.FunctionInCommon as FIC


USAGE = """
Usage: MFEprimer [options]

    Example 1: MFEprimer -i input_file -d database_file

    Example 2: MFEprimer -i input_file -d database_file -o output_file -W 7 -e 1000

    Example 3: MFEprimer.py -i primers.fa -d ./database.seq -T F -e 10000 -W 11 -a 2 --ppc_cutoff=0.2 --seq=primer.mfe.fa -o primer.mfe

Options:

Running 'MFEprimer -h' for USAGE information

    -h  Print this USAGE information

These options are mandatory for MFEprimer and are set by the user 

    -i  Query File [File In], Needed
    -d  Database [File In], Needed

These options are specified and optional for MFEprimer

    -o  MFEprimer Output File [File Out]
	default will print the results to stdout
    -s  The cutoff value of PPC for the MFEprimer output [Float]
	default = 0 (A float value between 0 and 1)
	Example: -s 0.4 (Means that only amplicons which PPC value is higher 
	than 0.4 would be reported on the output file)
    -ppc_cutoff The same as -s, for easy to memory
    -b	Beginning of the expected amplicon size range in bp[Integer]
        default = 0
    -B  Ending of the expected amplicon size range in bp[Integer]
        default = 5000
    -T  Format the database use 'formatdb -d database_file -p F -o T' [T/F] 
	default = F (Choose T, if the database has not been formated already)
    -seq Output the amplicons sequences in Fasta format
    -fasta The same as -seq, Output the amplicons sequences in Fasta format

These options are for calculating the Tm and Gibbs free energy
    --monovalent (mM)    Concentration of monovalent cations, default is 50.0 mM 
    --divalent (mM)      Concentration of divalent cations, default is 1.5 mM
    --oligo (nM)         Annealing oligo concentration, default is 50 nM
    --dNTP (mM)          Concentration of dNTPs, default is 0.25 mM

Below are some BLAST arguments available in MFEprimer

    -e  Expectation value (E) [Real]
        default = 1000.0
    -F  Filter query sequence (DUST with blastn, SEG with others) [String]
        default = F 
    -G  Cost to open a gap (-1 invokes default behavior) [Integer]
        default = -1
    -E  Cost to extend a gap (-1 invokes default behavior) [Integer]
        default = -1
    -q  Penalty for a nucleotide mismatch (blastn only) [Integer]
        default = -3
    -r  Reward for a nucleotide match (blastn only) [Integer]
        default = 1
    -f  Threshold for extending hits 
        default = 0 
    -g  Perform gapped alignment (not available with tblastx) [T/F]
        default = T
    -Q  Query Genetic code to use [Integer]
        default = 1
    -a  Number of processors to use [Integer]
        default = 1
    -M  Matrix [String]
        default = BLOSUM62
    -W  Word size, default if zero (blastn 11, megablast 28, all others 3) [Integer]
        default = 0
"""

def checkOpt_doBlast():
    """
    Checking and validating parameters
    Do BLASTN
    Return the parameters and Blast report file
    """
    # Check opt
    if len(sys.argv) <= 1:
        print >> sys.stderr, USAGE
        exit(1)
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], \
                    "i:d:o:s:b:B:e:F:G:E:q:v:f:g:Q:a:M:W:T:h", ["--help", "monovalent=", "divalent=", "oligo=", "dNTP=", "ppc_cutoff=", "seq=", "fasta="])
        #Essential arguments
        arg_i = '';
        arg_d = '';
        #Set default values for the optional arguments
        ppc_cutoff = 0
        range_start = 0
        range_stop = 5000
        format_wish = 'F'
        output = 'STDOUT'; # For MFEprimer output NOT for blastall
        # Default print the results to STDOUT
        seq_fasta_out = ''

        # For calculating Tm and Gibbs free energy
        monovalent_conc = 50.0
        divalent_conc = 1.5
        oligo_conc = 50.0
        dNTP_conc = 0.25

        # default value in BLASTN
        arg_e = 1000
        arg_F = 'F'
        arg_G = -1
        arg_E = -1
        arg_q = -3
        arg_r = 1
        arg_f = 0
        arg_g = 'F'
        arg_Q = 1
        arg_a = 1
        arg_M = 'BLOSUM62'
        arg_W = 11

        for o, a in opts:
            #Print help infomation
            if o == "-h" or o == "--help":
                print USAGE
                sys.exit(1)

            #Checking some essential arguments
            if o == '--seq':
                seq_fasta_out = a

            if o == '--fasta':
                seq_fasta_out = a

            if o == '-i':
                arg_i = a

            if o == '-d':
                arg_d = a

            if o == '--monovalent':
                monovalent_conc = float(a)

            if o == '--divalent':
                divalent_conc = float(a) 

            if o == '--oligo':
                oligo_conc =  float(a)

            if o == '--dNTP':
                dNTP_conc =  float(a)
            #Checking MFEprimer only arguments 
            if o == '-o':
                output = a

            if o == '-s':
                ppc_cutoff = a
                ppc_cutoff = float(ppc_cutoff)
                ppc_cutoff = ppc_cutoff * 100

            if o == '--ppc_cutoff':
                ppc_cutoff = a
                ppc_cutoff = float(ppc_cutoff)
                ppc_cutoff = ppc_cutoff * 100

            if o == '-b':
                range_start = a
                range_start = int(range_start)

            if o == '-B':
                range_stop = a
                range_stop = int(range_stop)

            if o == '-T':
                format_wish = a

            #Checking BLAST usually arguments 
            if o == '-e':
                arg_e = a
                arg_e = int(arg_e)

            if o == '-F':
                arg_F = a

            if o == '-G':
                arg_G = a

            if o == '-E':
                arg_E = a

            if o == '-q':
                arg_q = a

            if o == '-r':
                arg_r = a

            if o == '-f':
                arg_f = a


            if o == '-g':
                arg_g = a

            if o == '-Q':
                arg_Q = a

            if o == '-a':
                arg_a = a

            if o == '-M':
                arg_M = a

            if o == '-W':
                arg_W = a

        if not arg_i:
            print >> sys.stderr, "The primer file in FASTA format is required!"
            print >> sys.stderr, USAGE
            exit(1)

        if not arg_d:
            print >> sys.stderr, "Database file is required!"
            print >> sys.stderr, USAGE
            exit(1)

    except getopt.GetoptError:
        print >> sys.stderr, USAGE
        exit(1)

    # Format database
    if format_wish == 'T':
        format_cmd = '%s/formatdb -i %s -p F -o T' % (blast_home, arg_d)
        subprocess.Popen(format_cmd, shell=True).wait()

    # Do blast

    blast_file = '%s%s.%s.%s.MFEprimer.tmp' % (os.curdir, os.sep, time.time(), FIC.random_string(20))
    blast_cmd = '%s/blastall -p blastn -i %s -d %s -o %s -e %s -F %s -G %s -E %s -q %s -r %s -f %s -g %s -Q %s -a %s -M %s -W %s -b 1000000 -v 0' \
            % (blast_home, arg_i, arg_d, blast_file, arg_e, arg_F, arg_G, arg_E, arg_q, arg_r, \
               arg_f, arg_g, arg_Q, arg_a, arg_M, arg_W)

    # checking command 
    if not os.path.exists(blast_cmd.split()[0]):
        print >> sys.stderr, "blastall does not exist at %s" % blast_cmd

    subprocess.Popen(blast_cmd, shell=True).wait()
    # This will be printed in the end of the MFEprimer output file
    page_end = """
#####################################
#   Arguments used in MFEprimer   #
#####################################
    """
    page_end = page_end + '\nAmplicons size range from: %s bp to: %s bp\n' % (range_start, range_stop)
    page_end = page_end + 'PPC cutoff value: %s%s\n\n' % (ppc_cutoff, '%')
    #page_end = page_end + 'Main arguments used in BLASTN:\n\n'
    page_end = page_end + 'Expectation value: %s\n' % arg_e
    page_end = page_end + 'Word size: %s\n\n' % arg_W
    #page_end = page_end + 'Arguments used for calculating Tm and Gibbs free energy\n\n'
    page_end = page_end + 'Concentration of monovalent cations (usually KCl): %s mM\n' % monovalent_conc
    page_end = page_end + 'Concentration of divalent cations (usually MgCl2): %s mM\n' % divalent_conc
    page_end = page_end + 'Annealing oligo concentration: %s nM\n' % oligo_conc
    page_end = page_end + 'Concentration of dNTPs: %s mM\n' % dNTP_conc


    return (arg_i, output, blast_file, arg_d, ppc_cutoff, \
            range_start, range_stop, page_end, monovalent_conc, divalent_conc, oligo_conc, dNTP_conc, arg_W, seq_fasta_out)

def parseFasta(arg_i):
    """
    Parsing query file
    Return a dict containing the primers infomation
    """
    primer_info = {}
    handle = open(arg_i, 'r')
    sn = 0
    for seq_record in LFP.parse(handle):
        sn = sn + 1
        primer_info[seq_record.id] = {
            'seq' : seq_record.sequence, 
            'len' : seq_record.length, 
            'sn'  : sn,
            }
    handle.close()

    return primer_info

def reverseString(s):
    """Reverses a string given to it."""
    return s[::-1]


def parseBlast(blast_file, primer_info, word_size):
    """
    Parsing Blast output file
    Return plus_info dict which contains condidate forward primers;
    minus_info dict which contains condidate reverse primers;
    title_info dict for MFEprimer output 
    """
    #Parse BLAST output
    fh = open(blast_file, 'r')
    result = LBP.parse(fh)
    fh.close()
    #os.remove(blast_file) # Clean the temp file

    plus_info = []
    minus_info = []
    title_info = {}

    for record in result.records:
        query_id = record.query_id
        if len(title_info) < 1:
            title_info = {
                'version'     : record.version, 
                'reference'   : record.reference, 
                'database'    : record.database, 
                'db_sequences' : record.db_sequences, 
                'db_letters'   : record.db_letters, 
                }
        binding_num = 0
        for alignment in record.hits:    
            hit_id = alignment.hit_id
            binding_num = binding_num + len(alignment.hsps)
            for hsp in alignment.hsps: #Hsps (may be more)
                query_length = record.query_letters
                sbjct_strand = hsp.sbjct_strand

                # Collecting
                if sbjct_strand == 'Minus':
                    hsp_info = filterHsp(hsp, query_id, \
                                         hit_id, query_length, word_size)
                    if hsp_info:
                        minus_info.append(hsp_info)
                else:
                    hsp_info = filterHsp(hsp, query_id, \
                                         hit_id, query_length, word_size)
                    if hsp_info:
                        plus_info.append(hsp_info)                    

        primer_info[query_id]['binding_num'] = binding_num                

    return (title_info, plus_info, minus_info, primer_info)



def filterGap(qseq, aseq, sseq, sbjct_end):
    """
    Handle the HSP which alignment have gaps made by insert space
    """
    # If alignment in HSP have gaps made by insert space etc,
    # This would not thappen in matching between primer and tmplate for
    # real PCR reaction, so, we splited the alignment based on '-' and 
    # removed the part near 5'

    qpos = qseq.find('-')
    spos = sseq.find('-')
    if qpos > -1 and spos > -1: # Actually pos value should greater than 2
        if qpos <= spos:
            pos = qpos
            aseq = aseq[:pos] 
            qseq = qseq[:pos] + qseq[(pos + 1) :].lower()
            sseq = sseq[:pos]
            sbjct_end = sbjct_end - len(sseq)
            return (qseq, aseq, sseq, sbjct_end)
        else:
            pos = spos
            qseq = qseq[:pos] + qseq[pos:].lower()
            sseq = sseq[:pos]
            aseq = aseq[:pos]
            sbjct_end = sbjct_end - len(sseq)
            return (qseq, aseq, sseq, sbjct_end)
    elif qpos > -1 and spos == -1:
        pos = qpos
        aseq = aseq[:pos]
        qseq = qseq[:pos] + qseq[(pos + 1) :].lower()
        sseq = sseq[:pos]
        sbjct_end = sbjct_end - len(sseq)
        return (qseq, aseq, sseq, sbjct_end)
    elif spos > -1 and qpos == -1:
        pos = spos
        qseq = qseq[:pos] + qseq[pos:].lower()
        sseq = sseq[:pos]
        aseq = aseq[:pos]
        sbjct_end = sbjct_end - len(sseq)
        return (qseq, aseq, sseq, sbjct_end)
    else:
        return (qseq, aseq, sseq, sbjct_end)

def filterHsp(hsp, query_id, hit_id, query_length, word_size):
    """ Filter the HSPs which have mismatches in 3'"""

    sbjct_strand = hsp.sbjct_strand
    query_end   = hsp.query_end

    if sbjct_strand == 'Minus':
        aseq = reverseString(hsp.aseq)

        if query_length == query_end and len(aseq) >= int(word_size) and aseq[:int(word_size)].count(' ') < 1:
            sbjct_begin = hsp.sbjct_end # exchange
            sbjct_end   = hsp.sbjct_begin # exchange

            qseq = LS.seq(hsp.qseq).rev_com
            sseq = LS.seq(hsp.sseq).rev_com
            qseq, aseq, sseq, sbjct_end = \
                    filterGap(qseq, aseq, sseq, sbjct_end)
            hsp_info = { 
                    'query_id' : query_id,
                    'hit_id' : hit_id,
                    'query_length'   : query_length, 
                    'query_begin' : hsp.query_begin, 
                    'query_end'   : query_end, 
                    'sbjct_begin' : sbjct_begin, 
                    'sbjct_end'   : sbjct_end, 
                    'qseq'    : qseq, 
                    'aseq'    : aseq, 
                    'sseq'    : sseq, 
                    'score' : hsp.score,
                    }
            return hsp_info
        else:
            return 0
    else:
        aseq = hsp.aseq
        if query_length == query_end and len(aseq) >= int(word_size) and aseq[-int(word_size):].count(' ') < 1:
            sbjct_begin = hsp.sbjct_begin
            sbjct_end   = hsp.sbjct_end
            qseq = hsp.qseq
            sseq = hsp.sseq 
            qseq, aseq, sseq, sbjct_end = \
                    filterGap(qseq, aseq, sseq, sbjct_end)
            hsp_info = { 
                    'query_id' : query_id,
                    'hit_id' : hit_id,
                    'query_length'   : query_length, 
                    'query_begin' : hsp.query_begin, 
                    'query_end'   : query_end, 
                    'sbjct_begin' : sbjct_begin, 
                    'sbjct_end'   : sbjct_end, 
                    'qseq'    : qseq, 
                    'aseq'    : aseq, 
                    'sseq'    : sseq, 
                    'score' : hsp.score,
                    }            
            return hsp_info
        else:
            return 0



def analysis(primer_info, plus_info, minus_info, \
             range_start, range_stop, ppc_cutoff, \
             db_file, monovalent_conc, divalent_conc, \
             oligo_conc, dNTP_conc):
    """ 
    Compared one by one to find out the proper forward and reverse primers
    """
    result = []
    for p in plus_info:  # phsp stand for plus phsp
        for m in minus_info:
            p_hit_id = p['hit_id']
            m_hit_id = m['hit_id']
            if p_hit_id == m_hit_id:
                amplicon = compare(p_hit_id, p, m, range_start, \
                                   range_stop, ppc_cutoff, \
                                   primer_info, db_file, \
                                   monovalent_conc, divalent_conc, \
                                   oligo_conc, dNTP_conc)

                if amplicon:
                    result.append(amplicon)

    return result

def compare(hit_id, p, m, range_start, range_stop, \
            ppc_cutoff, primer_info, database, \
            monovalent_conc, divalent_conc, \
            oligo_conc, dNTP_conc):
    """
    Compared according to real PCR reaction
    """
    p_query_id = p['query_id']
    p_query_length = p['query_length']
    p_query_begin = p['query_begin']
    p_query_end = p['query_end']
    p_sbjct_begin = p['sbjct_begin']
    p_sbjct_end = p['sbjct_end']
    p_qseq = p['qseq']
    p_aseq = p['aseq']
    p_sseq = p['sseq']
    p_score = p['score']
    
    m_query_id = m['query_id']
    m_query_length = m['query_length']
    m_query_begin = m['query_begin']
    m_query_end = m['query_end']
    m_sbjct_begin = m['sbjct_begin']
    m_sbjct_end = m['sbjct_end']
    m_qseq = m['qseq']
    m_aseq = m['aseq']
    m_sseq = m['sseq']
    m_score = m['score']
    
    amplicon = 0
    if m_sbjct_begin >= p_sbjct_end:
        length = m_sbjct_begin - p_sbjct_end + \
                p_query_length + m_query_length - 1

        if length >= range_start and length <= range_stop:
            ppc = cal_ppc(p_aseq, m_aseq, p_query_length, m_query_length)

            if ppc > ppc_cutoff:
                amplicon = formatShow(primer_info, hit_id, \
                                      p_query_id, m_query_id, \
                                      p_sbjct_end, m_sbjct_begin, \
                                      p_qseq, p_sseq, m_qseq, \
                                      m_sseq, p_query_begin,\
                                      m_query_begin, p_aseq, m_aseq, \
                                      length, ppc,\
                                      database, \
                                      p_score, m_score, \
                                      monovalent_conc, divalent_conc, \
                                      oligo_conc, dNTP_conc)

    
    return amplicon            
    
def formatShow(primer_info, hit, fid, rid, f3, r3, fqseq, fhseq, rqseq, \
               rhseq, fstart, rstart, \
               fmatch, rmatch, length, ppc, \
               database, 
               f_score, r_score, \
               monovalent_conc, divalent_conc, \
               oligo_conc, dNTP_conc):
    """
    Format the results for friendly output
    """
    hit_tmp = re.sub('\.\d+\|', '|', hit)
    hit_tmp = re.sub('\|', '\\|', hit_tmp)
    
    # Extract sequence
    if r3 > f3 + 2:
        cmd = '%s/fastacmd -s %s -d %s -L %d,%d' % (blast_home, hit_tmp, database, (f3 + 1), (r3 - 1))
    
        (outtext, outerror) = subprocess.Popen(cmd, shell=True, 
                                               stdout=subprocess.PIPE).communicate()
    
        if outerror:
            print >> sys.stderr, "fastacmd running error!"
            sys.exit(1)
    
        amplicon_def = string.split(outtext, '\n')[0]
        amplicon_def = ' '.join(string.split(amplicon_def, ' ')[1:])
        amplicon_seq = ''.join(string.split(outtext, '\n')[1:])
    
    else:
        cmd = '%s/fastacmd -s %s -d %s -L %d,%d' % (blast_home, hit_tmp, database, f3, r3)
    
    
        (outtext, outerror) = subprocess.Popen(cmd, shell=True, \
                                               stdout=subprocess.PIPE).communicate()
    
        if outerror:
            print >> sys.stderr, "%fastacmd running error!"
            sys.exit(1)
    
        amplicon_def = string.split(outtext, '\n')[0]
        amplicon_def = ' '.join(string.split(amplicon_def, ' ')[1:])
        amplicon_seq = ''
    
    amplicon_seq = string.lower(amplicon_seq)
    if re.match('No definition line found', amplicon_def):
        amplicon_def = ''
    else:
        amplicon_def = re.sub('>', '', amplicon_def)
    
    ppc_tmp = '%.1f' % ppc
    if fid == rid:
        ppc_tmp = '-' + ppc_tmp
    
    #line_0 = '  PPC = %s%s\n' % (ppc_tmp, '%')
    
    f_tm = TD.calTm(fqseq, fhseq, monovalent_conc, divalent_conc, \
                 oligo_conc, dNTP_conc)
    if f_tm < 0:
    	f_tm_tmp = '-'
    else:
    	f_tm_tmp = '%.1f' % f_tm

    if len(fqseq) < 5:
        f_deltaG = TD.calDeltaG(fqseq, fhseq, monovalent_conc, divalent_conc, dNTP_conc)
        f_deltaG_3 = f_deltaG
    else:
        f_deltaG_3 = TD.calDeltaG(fqseq[-5:], fhseq[-5:], monovalent_conc, divalent_conc, dNTP_conc)
        f_deltaG = TD.calDeltaG(fqseq, fhseq, monovalent_conc, divalent_conc, dNTP_conc) 
    
    line_1 = '  FP: 3\'ΔG = %.1f%s, Tm = %s %s\n' \
            % (f_deltaG_3, '(kcal/mol)', f_tm_tmp, '(°C)')

    r_tm = TD.calTm(rqseq, rhseq, monovalent_conc, divalent_conc, \
                 oligo_conc, dNTP_conc)
    if r_tm < 0:
    	r_tm_tmp = '-'
    else:
    	r_tm_tmp = '%.1f' % r_tm

    if len(rqseq) < 5:
        r_deltaG = TD.calDeltaG(rqseq, rhseq, monovalent_conc, divalent_conc, dNTP_conc)
        r_deltaG_3 = f_deltaG
    else:
	r_deltaG_3 = TD.calDeltaG(rqseq[:5], rhseq[:5], monovalent_conc, divalent_conc, dNTP_conc)
        r_deltaG = TD.calDeltaG(rqseq, rhseq, monovalent_conc, divalent_conc, dNTP_conc) 
    
    line_1 = line_1 + '  RP: 3\'ΔG = %.1f%s, Tm = %s %s\n' \
            % (r_deltaG_3, '(kcal/mol)', r_tm_tmp, '(°C)')
    
    fhead = primer_info[fid]['seq'][0 : (fstart - 1)]
    fhead = fhead.lower()
    rtail = primer_info[rid]['seq'][0 : (rstart - 1)]
    rtail = rtail.lower()
    rtail = LS.seq(rtail).rev_com # Thanks Prof. Zhang for find this bug
    
    line_3 = '  Binding sites: %s(%s/%s) ... %s(%s/%s)\n\n' % \
            ((f3 - len(fhseq) + 1 - fhseq.count(' ')), \
             (len(fmatch) - fmatch.count(' ')), \
             primer_info[fid]['len'], \
             (r3 + len(rhseq) - 1), \
             (len(rmatch) - rmatch.count(' ')), \
             primer_info[rid]['len'])
    
    seq4gc = fhead + fqseq + amplicon_seq + rqseq + rtail
    GC_content = cal_GC(seq4gc, length)

    
    line_0 = '  PPC = %s%s, Size = %s bp, GC content = %.1f%s\n' % (ppc_tmp, '%', length, GC_content, '%')
    
    tmp_line = '>>>' + fid
    line_5 = tmp_line.ljust(0) + '\n'
    
    if fstart > 1:
        line_6 = '%s%s%s%s%s%s\n' % (3 * ' ', '1', (fstart - 2) * ' ' , fstart , \
                                     (len(fqseq) - len(str(fstart)) - 1) * ' ',\
                                      primer_info[fid]['len'])
    else:
        line_6 = '%s%s%s%s\n' % (3 * ' ', '1', \
                                 (len(fqseq) - 2) * ' ', primer_info[fid]['len'])
    
    line_7 = '5\' ' + fhead + fqseq + ' 3\'\n'
    len_1 = len(fhead) + len(fqseq) + len(rqseq) + len(rtail) + 9

    if len_1 <= 74:
        len_2 = 80
    else:
        len_2 = len_1 + 6

    middle_length = int((len_2 - len_1) / 2)
    
    line_9 = len(fhead) * ' ' + '5\' ' + fhseq + amplicon_seq[0 : middle_length] \
                + '...' + amplicon_seq[-middle_length : ] + rhseq + ' 3\''
    line_8 = '%s%s%s%s\n' % ((len(fhead) + 3) * ' ', 
                             fmatch, \
                             (len(amplicon_seq[0:middle_length]) \
                              + len(amplicon_seq[-middle_length:]) \
                              + 3 + len(rhseq) - 1) * ' ', \
                              (r3 + len(rhseq) - 1 - rhseq.count(' ')))
    line_10 = '%s%s%s%s\n' % ((3 + len(fhead)) * ' ', \
                              (f3 - len(fhseq) + 1 + fhseq.count(' ')), \
                              (len(fhseq) + (len(amplicon_seq[0:middle_length]) \
                                + len(amplicon_seq[-middle_length:])) \
                                + 3 - len(str(f3 - len(fhseq) + 1))) * ' ', \
                                rmatch) 
    line_11 = (len(line_9) - len(rhseq) - 6) * ' ' + '3\' ' \
                + rqseq + rtail + ' 5\'' + '\n'
    
    if rstart > 1:
        line_12 = '%s%s%s%s%s%s\n' % ((len(line_9) - 3 - len(rhseq)) * ' ', \
        primer_info[rid]['len'], \
        (len(rhseq) - len(str(primer_info[rid]['len'])) - 1) * ' ', rstart, \
        (len(rtail) - len(str(rstart))) * ' ', '1')
    else:
        line_12 = '%s%s%s%s\n' % ((len(line_9) - 3 - len(rhseq)) * ' ', \
        primer_info[rid]['len'], \
        (len(rqseq) - 1 - len(str(primer_info[rid]['len']))) * ' ', '1')
    
    tmp_line = rid + '<<<'
    line_13 = tmp_line.rjust(80) + '\n\n'
    
    line = line_0 + line_1 + line_3 + line_5 + line_6 \
            + line_7 + line_8 + line_9 + '\n' + line_10 + line_11 \
            + line_12 + line_13
    
    
    seq = FIC.print_seq(seq4gc, width = 76) + '\n\n\n'
    
    amplicon = {
        'geneid' : hit, 
        'fid' : fid, 
        'rid' : rid, 
        'size': length, 
        'ppc' : ppc, 
        'desc'   : amplicon_def, 
        'line'   : line, 
        'seq'    : seq,

        'fp_deltaG' : f_deltaG,
        'fp_3_deltaG' : f_deltaG_3,
        'fp_tm' : f_tm,

        'rp_deltaG' : r_deltaG,
        'rp_3_deltaG' : r_deltaG_3,
        'rp_tm' : r_tm,
        }
    
    return amplicon

    

def display_title(primer_info, title_info, output):
    """
    Display the title part of MFEprimer output
    """

    line = ' '.join(sys.argv) + '\n\n'

    line = line + 'MFEprimer %s [%s]\n%s\n\n' % (__version__, __date__, __reference__)

    #line = line + '%s %s [%s]\n\n' % (title_info['application'], title_info['version'], title_info['date'])
    line = line + title_info['version'] + '\n'

    tmp_line = '%s\n\n' % title_info['reference']
    tmp_line = re.sub('~', '', tmp_line, 1)
    tmp_line = re.sub('~', ' ', tmp_line)
    line = line + textwrap.fill(tmp_line, width = 74)

    tmp_line = ''
    seq_num = 0
    letters = 0
    for primer_id in primer_info.keys():
        tmp_line = tmp_line + primer_id + '; '
        seq_num = seq_num + 1
        letters = letters + primer_info[primer_id]['len']
    tmp_line = re.sub(', $', '', tmp_line)

    line = line + '\n\n\nQuery = %s\n' % tmp_line
    line = line + '        (%s primer sequences; %s letters)\n\n' % (seq_num, letters)

    line = line + 'Database: %s\n' % title_info['database']
    line = line + '           %s sequences; %s total letters\n\n' % \
            (title_info['db_sequences'], title_info['db_letters'])
    line = line + 'Searching..................................................done\n\n\n'

    if output == 'STDOUT':
        print line
    else:
        f = file(output, 'w') 
        f.write(line)
        f.close

def cal_ppc(f_match, r_match, f_size, r_size):
    """
    Calculate the PPC value
    """
    # Many thanks Prof. Zhang for advices of PPC calculating
    f_match = len(f_match) - f_match.count(' ')
    r_match = len(r_match) - r_match.count(' ')
    ave = (f_match + r_match) / 2
    stdev = math.sqrt((f_match - ave)**2 + (r_match - ave)**2)
    cv = 1 - stdev/ave
    ppc = f_match / f_size * r_match / r_size * cv * 100

    return ppc



def display(bin, output, seq_fasta_out):
    """
    Display the main part of MFEprimer output
    """
    num = len(bin)
    line = 'Distribution of %s MFEprimer hits on the query primers\n\n' % num
    line = line + '[Sort by primer pair coverage (PPC)]\n'
    line = line + 'FP'.rjust(89)+ 'RP'.rjust(6) + 'FP'.rjust(8)+ 'RP'.rjust(6) + '\n'
    line = line + 'Size'.rjust(73) + 'PPC'.rjust(9) + '3\'ΔG'.rjust(9) + '3\'ΔG'.rjust(7) + 'Tm'.rjust(7) + 'Tm'.rjust(6) + '\n'
    line = line + 'Sequences producing potential PCR products:'.ljust(63) + '(bp)'.rjust(10) + '(%)'.rjust(9) + '(kcal/mol)'.rjust(14) + '(°C)'.rjust(9) + '(°C)'.rjust(7) + '\n\n'


    # Sorting the results first by ppc and then by size
    list4sort = []
    for th in range(num):
        tmp = (th, bin[th]['ppc'], bin[th]['size'])
        list4sort.append(tmp)

    list4sort.sort(lambda x, y: (
        cmp(y[1], x[1])
        or
        cmp(x[2], y[2])))

    ref = {}
    sn = 0
    for item in list4sort:
        sn = sn + 1
        th = item[0]
        ref[th] = sn
        desc = bin[th]['desc']
        geneid = bin[th]['geneid']
        desc = '%s: %s %s' % (sn, geneid, desc)
        if len(desc) > 44:
            desc = desc[:44] + '...'

        size = bin[th]['size']
        fid = bin[th]['fid']
        rid = bin[th]['rid']
        ppc = bin[th]['ppc']
        ppc = '%.1f' % ppc
	fp_3_deltaG = '%.1f' % bin[th]['fp_3_deltaG']
	rp_3_deltaG = '%.1f' % bin[th]['rp_3_deltaG']
	fp_tm = '%.1f' % bin[th]['fp_tm']
	rp_tm = '%.1f' % bin[th]['rp_tm']

        if fid == rid:
            ppc = '-' + ppc

	line = line + desc.ljust(69) + str(size).ljust(7) + ' '*2 + str(ppc).ljust(6) + ' '*2 + str(fp_3_deltaG).ljust(5) + ' '*2 + str(rp_3_deltaG).ljust(5) + ' '*2 + str(fp_tm).ljust(4) + ' '*2 + str(rp_tm).ljust(4) + '\n'


    line = line + '\n\n'

    line = line + '[Sort by amplicon size]\n'
    line = line + 'FP'.rjust(89)+ 'RP'.rjust(6) + 'FP'.rjust(8)+ 'RP'.rjust(6) + '\n'
    line = line + 'Size'.rjust(73) + 'PPC'.rjust(9) + '3\'ΔG'.rjust(9) + '3\'ΔG'.rjust(7) + 'Tm'.rjust(7) + 'Tm'.rjust(6) + '\n'
    line = line + 'Sequences producing potential PCR products:'.ljust(63) + '(bp)'.rjust(10) + '(%)'.rjust(9) + '(kcal/mol)'.rjust(14) + '(°C)'.rjust(9) + '(°C)'.rjust(7) + '\n\n'
    # Sorting the results first by size and than by ppc
    list4sort_2 = []
    for th in range(num):
        tmp = (th, bin[th]['size'], bin[th]['ppc'])
        list4sort_2.append(tmp)
    list4sort_2.sort(lambda x, y: (
        cmp(x[1], y[1])
        or
        cmp(y[2], x[2])))

    for item in list4sort_2:
        th = item[0]
        sn = ref[th]
        desc = bin[th]['desc']
        geneid = bin[th]['geneid']
        desc = '%s: %s %s' % (sn, geneid, desc)
        if len(desc) > 44:
            desc = desc[:44] + '...'

        size = bin[th]['size']
        fid = bin[th]['fid']
        rid = bin[th]['rid']
        ppc = bin[th]['ppc']
        ppc = '%.1f' % ppc
	fp_3_deltaG = '%.1f' % bin[th]['fp_3_deltaG']
	rp_3_deltaG = '%.1f' % bin[th]['rp_3_deltaG']
	fp_tm = '%.1f' % bin[th]['fp_tm']
	rp_tm = '%.1f' % bin[th]['rp_tm']

        if fid == rid:
            ppc = '-' + ppc

	line = line + desc.ljust(69) + str(size).ljust(7) + ' '*2 + str(ppc).ljust(6) + ' '*2 + str(fp_3_deltaG).ljust(5) + ' '*2 + str(rp_3_deltaG).ljust(5) + ' '*2 + str(fp_tm).ljust(4) + ' '*2 + str(rp_tm).ljust(4) + '\n'

    line = line + '\n\n'

    line = line + "Details for the primers binding to the DNA template\n"
    line = line + "[Sort by primer pair coverage (PPC)]\n\n"
    seq_out_list = [] # For out put the amplicon sequence in Fasta sequence
    sn = 0
    for item in list4sort:
        sn = sn + 1
        th = item[0]
        geneid = bin[th]['geneid']
        desc = bin[th]['desc']
        fid = bin[th]['fid']
        rid = bin[th]['rid']
        seq = bin[th]['seq']
        line_part = bin[th]['line']
        line = line + '%s: %s + %s ==> %s %s\n\n' % (sn, fid, rid, geneid, desc)
        line = line + line_part
        tmp_line = '>%s_%s_%s_%s %s' % (sn, fid, rid, geneid, desc)
        if len(tmp_line) > 77:
            tmp_line = tmp_line[:77] + '...\n'
        else:
            tmp_line = tmp_line + '\n'

        seq_out_list.append(tmp_line)
        seq_out_list.append(seq.strip())
        seq_out_list.append('\n')
        line = line + tmp_line + seq

    if output == 'STDOUT':
        print line
    else:
        f = file(output, 'a') 
        f.write(line)
        f.close

    if seq_fasta_out:
        fh = file(seq_fasta_out, 'w')
        for line in seq_out_list:
            fh.write(line)
        fh.close

def cal_GC(seq, length):
    """
    Calculate the GC content of a sequence
    """
    seq = seq.lower()
    c_num = seq.count('c')
    g_num = seq.count('g')
    n_num = seq.count('n')

    GC_value = (c_num + g_num + n_num * 0.5) / length * 100
    return GC_value

def display_end(output, page_end):
    """
    Display the end of the MFEprimer output
    """
    if output == 'STDOUT':
        print page_end
    else:
        f = file(output, 'a') 
        f.write(page_end)
        f.close

def main():

    # Checking parameters and then Do Blast
    (arg_i, output, blast_file, arg_d, ppc_cutoff, range_start, \
     range_stop, page_end, monovalent_conc, divalent_conc, oligo_conc, dNTP_conc, word_size, seq_fasta_out) \
    = checkOpt_doBlast()
    
    # Parsing query file
    primer_info = parseFasta(arg_i)
   
    # Parsing Blast reprot file 
    (title_info, plus_info, minus_info, primer_info) = parseBlast(blast_file, primer_info, word_size)
    
    # Analyzing and find out the proper forward and reverse primers
    result = analysis(primer_info, plus_info, minus_info, \
                      range_start, range_stop, ppc_cutoff, \
                      arg_d, monovalent_conc, divalent_conc, \
                      oligo_conc, dNTP_conc)

    # Displaying
    display_title(primer_info, title_info, output)
    display(result, output, seq_fasta_out)
    display_end(output, page_end)
    os.remove(blast_file)

if __name__ == '__main__':
    main()
