#!/usr/bin/python3

'''
Title: sprot_extractID.py
Date: 2018-01-11
Author: Adam Nunn

Description:
    This program works alongside "sprot_fetchID.py" to obtain specific
    sequences from the uniprot_sprot.pep database, eg. all sequences
    from Fungi or Diptera. The program takes two input files: 1) the
    list of sequences to be extracted from the database (INPUT_IDS_FILE)
    and 2) the uniprot_sprot database in protein fasta format
    (INPUT_FASTA_FILE). This program will simply extract the sequences
    specified in INPUT_IDS_FILE from the database file INPUT_FASTA_FILE
    and output to a new protein fasta file (OUTPUT_FASTA_FILE).

Usage:
    ./sprot_extractID.py INPUT_IDS_FILE INPUT_FASTA_FILE OUTPUT_FASTA_FILE
eg. ./sprot_extractID.py Fungi.txt uniprot_sprot.pep Fungi.uniprot_sprot.pep

'''

import sys

idList = []
flag = False

with open(sys.argv[1], 'r') as INPUT_IDS_FILE, open(sys.argv[2], 'r') as INPUT_FASTA_FILE, open(sys.argv[3], 'w') as OUTPUT_FASTA_FILE:
 for line in INPUT_IDS_FILE: # parse through the sequences IDs
  line = line.rstrip()
  idList.append(line)
 for sequence in INPUT_FASTA_FILE: # parse through the uniprot_sprot.pep database
  sequence = sequence.rstrip()   
  if sequence.startswith(">"):
   sequence = sequence[1:]
   ID = sequence.split(" ")
   if ID[0] in idList:
    flag = True
    print(">{}".format(sequence), file=OUTPUT_FASTA_FILE) # print found sequence ID to output
   else: flag = False
  else:
   if flag == True:
    print(sequence, file=OUTPUT_FASTA_FILE) # print corresponding sequence to output
   else: continue
