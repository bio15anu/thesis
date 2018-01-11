#!/usr/bin/python3

'''
Title: sprot_fetchID.py
Date: 2018-01-11
Author: Adam Nunn

Description:
    This program works alongside "sprot_extractID.py" to obtain specific
    sequences from the uniprot_sprot.pep database, eg. all sequences
    from Fungi or Diptera. The program takes a single input file: 1) the
    uniprot_sprot.dat file containing metadata for each sequence in the
    uniprot_sprot.pep database (INPUT_DAT_FILE). This program will
    obtain sequence IDs for all sequences from a user-specified 
    classification (ORGANISM), eg. Fungi NOTE: the naming convention for
    this user command must match the naming convention in the
    uniprot_sprot.dat file (without the ";" character). The obtained
    sequence IDs will be output to a file in the current directory
    (OUTPUT_FASTA_FILE), to be used in the follow-up program.

Usage:
    ./sprot_fetchID.py INPUT_DAT_FILE OUTPUT_IDS_FILE ORGANISM
eg. ./sprot_fetchID.py uniprot_sprot.dat Fungi.txt Fungi

'''

import sys

ORGANISM = sys.argv[3]
idList = []
flag = False

with open(sys.argv[1], 'r') as INPUT_DAT_FILE, open(sys.argv[2], 'w') as OUTPUT_IDS_FILE:
 for line in INPUT_DAT_FILE: # parse through the uniprot_sprot.dat file
  line = line.rstrip()
  if line.startswith("ID") and (flag == False): # look for IDs
   line = line.split(" ")
   ID = line[3]
   flag = True
  elif line.startswith("OC") and (flag == True): # found ID, now look for ORGANISM
   line = line.split(" ")
   if (ORGANISM + ";") in line:
    print(ID, file=OUTPUT_IDS_FILE) # print sequence ID to output file
    flag = False
  elif line.startswith("//"): flag = False # reset flag, no ORGANISM data found
