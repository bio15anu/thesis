#!/usr/bin/python3

'''
Title: secretome_parseSSP.py
Date: 2018-01-11
Author: Adam Nunn

Description:
    This program takes a single input file: 1) a protein fasta file 
    containing seqeunces thought to be potential SSPs (INPUT_FASTA), eg.
    annotated with SignalP but not with TMHMM. The program will filter
    these sequences according to the specified proportion of cysteine
    residues (PROP), and the specified maximum length (LENG) to give an
    OUTPUT fasta file "NAME.PROP.LENGaa.fasta" in the current directory.

Usage:
    ./secretome_parseSSP.py INPUT_FASTA PROP LENG NAME
eg. ./secretome_parseSSP.py unfiltered.secretome.fasta 0.03 300 SSP

'''

import sys

PROP = sys.argv[2]
LENG = sys.argv[3]
NAME = sys.argv[4]

with open(sys.argv[1], 'r') as INPUT_FASTA, open(NAME + "." + PROP + "." + LENG + "aa.fasta", 'w') as OUTPUT:
 for line in INPUT_FASTA: # parse the input protein fasta file
  line = line.rstrip()
  if line.startswith(">"):
   seqID = line.rstrip()
  else:
   line = line.upper() # convert seqeunce to upper case   
   length = int(len(line)) # calculate length of sequence
   cyst = line.count("C") # calculate number of cysteine residues in sequence
   prop = int(cyst)/int(length) # calculate proportion of sequence made up of cysteine residues
   if (length <= int(sys.argv[3])) and (prop >= float(sys.argv[2])): # filter sequences and print to OUTPUT accordingly
    print(seqID, file=OUTPUT)
    print(line, file=OUTPUT)
