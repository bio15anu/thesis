#!/usr/bin/python3

'''
Title: alignments_filter_domains_from_lengths.py
Date: 2018-01-11
Author: Adam Nunn

Description:
    This program takes two input files: 1) tab-delimited file containing on each line a
    sequence ID and relative start/end positions for a pfam domain contained within the
    sequence (INPUT_LENGTHS_FILE) and 2) a fasta file containing the relevant sequence IDs
    and protein sequence on the following line (INPUT_FASTA_FILE). The output is a fasta
    file containing only the sequence IDs specified in the INPUT_LENGTHS_FILE and the
    extracted domain sequences on the following line (OUTPUT_FASTA_FILE).

Usage:
    ./alignments_filter_domains_from_lengths.py INPUT_LENGTHS_FILE INPUT_FASTA_FILE OUTPUT_FASTA_FILE
eg. ./alignments_filter_domains_from_lengths.py PF00089.lengths.txt Trinity_isoH_cdhit10.fasta.transdecoder.pep PF00089.fasta

'''

import sys
import re

outseqs = []

with open(sys.argv[1], 'r') as INPUT_LENGTHS_FILE, open(sys.argv[3], 'w') as OUTPUT_FASTA_FILE:
 for info in INPUT_LENGTHS_FILE: # parse the information in the INPUT_LENGTHS_FILE by line
  info = info.rstrip()
  info = info.split("\t")
  Found = False
  INPUT_FASTA_FILE = open(sys.argv[2], 'r') # open the INPUT_FASTA_FILE to find sequences
  for line in INPUT_FASTA_FILE:
   line = line.rstrip()
   if Found == False:
    matchObject = re.search(info[0] + ":", line)
    if matchObject: # sequence ID found in the INPUT_FASTA_FILE
     count = 1
     while (info[0] + "__" + str(count)) in outseqs: # check if sequence has already had domain extracted (some sequences contain domain copies)
      count += 1
     outseqs.append(info[0] + "__" + str(count)) # give new ID to sequence in case multiple domains are extracted from the same sequence
     print(">" + info[0] + "__" + str(count), file=OUTPUT_FASTA_FILE)
     Found = True
    else: continue # sequence ID not found, continue searching INPUT_FASTA_FILE
   else:
    start = int(info[1]) # define start position for domain
    stop = int(info[2]) # define end position for domain
    line = "".join(line[start:stop]) # extract domain from sequence
    print(line, file=OUTPUT_FASTA_FILE) # print domain sequence to OUTPUT_FASTA_FILE
    break
  INPUT_FASTA_FILE.close()
