#!/usr/bin/python3

'''
Title: sigP_mend_Transpep_for_sigP.py
Date: 2018-01-11
Author: Adam Nunn

Description:
    This program works alongside "sigP_mend_temp.py" to fix file
    formating issues when running SignalP on sequences with long
    sequence IDs. The program takes a single input file: 1) the
    original protein fasta file containing the sequences to be analysed
    and sequence IDs that are too long for SignalP (INPUT_FASTA_FILE).
    The program will simply rename the sequence IDs to a format that
    can be handled by SignalP, generating a new fasta file to be used
    with SignalP (OUTPUT_FASTA_FILE).

Usage:
    ./sigP_mend_temp.py INPUT_FASTA_FILE OUTPUT_FASTA_FILE
eg. ./sigP_mend_temp.py Trinity_isoH_cdhit10.fasta.transdecoder.pep signalp.transdecoder.pep

'''

import sys

count = 0

with open(sys.argv[1], 'r') as INPUT_FASTA_FILE, open(sys.argv[2], 'w') as OUTPUT_FASTA_FILE:
 for line in INPUT_FASTA_FILE:
  line = line.rstrip()
  if line.startswith(">"):
   print(">ID_{}".format(count), file=OUTPUT_FASTA_FILE)
   count += 1
  else: print(line, file=OUTPUT_FASTA_FILE)
