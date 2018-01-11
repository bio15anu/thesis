#!/usr/bin/python3

'''
Title: assembly_prefilter_with_Transpep.py
Date: 2018-01-11
Author: Adam Nunn

Description:
    This program takes two input files: 1) a fasta file containing predicted protein
    sequences, usually output from TransDecoder (INPUT_PEP_FASTA) and 2) a fasta file
    containing a nucleotide assembly to be filtered according to sequences present in
    the INPUT_PEP_FASTA file, usually output from Trinity (INPUT_ASSEMBLY_FASTA). The
    output is a fasta file of the nucelotide assembly that is filtered to contain only
    sequences also present in the INPUT_PEP_FASTA file (OUTPUT_ASSEMBLY_FASTA).

Usage:
    ./assembly_prefilter_with_Transpep.py INPUT_PEP_FASTA INPUT_ASSEMBLY_FASTA OUTPUT_ASSEMBLY_FASTA
eg. ./assembly_prefilter_with_Transpep.py Trinity_full.fasta.transdecoder.pep Trinity_full.fasta Trinity_trans.fasta

'''

import sys

transID = []
firstLine = True
flag = False

with open(sys.argv[1], 'r') as INPUT_PEP_FASTA, open(sys.argv[2], 'r') as INPUT_ASSEMBLY_FASTA, open(sys.argv[3], 'w') as OUTPUT_ASSEMBLY_FASTA:
 for info in INPUT_PEP_FASTA: # parse the INPUT_PEP_FASTA file
  info = info.rstrip()
  if info.startswith(">"):
   info = info.split(":")
   transID.append(info[2])
 transID = set(transID) # sequence IDs placed into a set
 
 for line in INPUT_ASSEMBLY_FASTA: # parse the INPUT_ASSEMBLY_FASTA file
  line = line.rstrip()
  if line.startswith(">"):
   listLine = line.split(" ")
   if (listLine[0][1:] in transID) and (firstLine == True): # find the first sequence in INPUT_ASSEMBLY_FASTA that is also in INPUT_PEP_FASTA
    print("{}".format(line), file=OUTPUT_ASSEMBLY_FASTA)
    flag = True
    firstLine = False
   elif (listLine[0][1:] in transID) and (firstLine == False): # find all but first sequence in INPUT_ASSEMBLY_FASTA that are also in INPUT_PEP_FASTA
    print("\n{}".format(line), file=OUTPUT_ASSEMBLY_FASTA)
    flag = True
   else: flag = False
  else:
   if flag == True:
    line = line.upper()
    print(line, end='', file=OUTPUT_ASSEMBLY_FASTA)
    
