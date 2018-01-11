#!/usr/bin/python3

'''
Title: GSEA_flip_rank_files.py
Date: 2018-01-11
Author: Adam Nunn

Description:
    This simple program takes a single input file: 1) a tab-delimited file
    containing sequences ranked according to a positive - negative scale
    (INPUT_RANK_FILE). The input rank file must be a temporary file, built
    from an existing rank file with positive scores at the top and negative
    scores at the bottom. The temporary file must contain the same data in
    reverse order. The program will simply then transform positive scores
    to negative and vice-versa in a new output rank file (OUTPUT_RANK_FILE).
    This is necessary for later use in GSEA analysis.

Usage:
    ./GSEA_flip_rank_files.py INPUT_RANK_FILE OUTPUT_RANK_FILE
eg. ./GSEA_flip_rank_files.py coni_versus_glen.tmp coni_versus_glen.rnk

'''

import sys

with open(sys.argv[1], 'r') as INPUT_RANK_FILE, open(sys.argv[1] + ".rnk", 'w') as OUTPUT_RANK_FILE:
 for line in INPUT_RANK_FILE: # parse the input rank file
  line = line.rstrip()
  line = line.split("\t")
  if line[1].startswith("-"):
   print(line[0] + "\t" + line[1][1:], file=rank_OUT) # remove "-" sign from negative scores
  else: print(line[0] + "\t-" + line[1], file=rank_OUT) # add "-" sign to positive scores
