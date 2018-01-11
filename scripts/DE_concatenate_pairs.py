#!/usr/bin/python3

'''
Title: DE_concatenate_pairs.py
Date: 2018-01-11
Author: Adam Nunn

Description:
    This program takes five input files: 1) a file containing sequence annotations,
    usually in either *.xls or *.gmt format (INPUT_ANNOTATIONS). The program will
    identify whether the input is in *.gmt format and if not it will assume *.xls
    format. The accepted *.xls file must be annotation report from Trinotate. Files
    2-5) output files from DESeq2 analysis, in each instance containing DE genes in
    one of four relevant pairwise comparisons NOTE: these files must be prefiltered
    to contain only DE genes upregulated in the specific condition of interest. The
    program will output a *.xls report file for the condition of interest,
    containing annotations and pval etc. for each pairwise comparison in one file.

Usage:
    ./DE_concatenate_pairs.py INPUT_ANNOTATIONS INPUT_PAIRWISE_DE1 INPUT_PAIRWISE_DE2 INPUT_PAIRWISE_DE3 INPUT_PAIRWISE_DE4 OUTPUT_FILE
eg. ./DE_concatenate_pairs.py Emuscae_annotation_report.xls coni_vs_glen coni_vs_symp coni_vs_para coni_vs_spor coni.xls

'''

import sys

# identify whether the input annotation file is in *.gmt format
if sys.argv[1][-3:] == "gmt": gmt_file = True
else: gmt_file = False

Id = []
Desc = []

pair1dict = {}
pair2dict = {}
pair3dict = {}
pair4dict = {}

with open(sys.argv[1], 'r') as INPUT_ANNOTATIONS, open(sys.argv[2], 'r') as INPUT_PAIRWISE_DE1, open(sys.argv[3], 'r') as INPUT_PAIRWISE_DE2, \
    open(sys.argv[4], 'r') as INPUT_PAIRWISE_DE3, open(sys.argv[5], 'r') as INPUT_PAIRWISE_DE4, open(sys.argv[6], 'w') as OUTPUT_FILE:
 pairTuple = (INPUT_PAIRWISE_DE1, INPUT_PAIRWISE_DE2, INPUT_PAIRWISE_DE3, INPUT_PAIRWISE_DE4)
 pairDictTuple = (pair1dict, pair2dict, pair3dict, pair4dict)
 for line in INPUT_ANNOTATIONS: # parse the INPUT_ANNOTATIONS FILE (*.xls or *.gmt)
  line = line.rstrip()
  line = line.split("\t")
  Id.append(line[0])
  if gmt_file == True: Desc.append(line[1]) # for GO terms etc. from a gmt file
  else: Desc.append(line[1:]) #for transcripts when using an xls file
 goDict = dict(zip(Id, Desc)) # create dictionary for sequence ID and annotations
 for i in range(4): # iterate through and parse each INPUT_PAIRWISE_DE file
  pair = pairTuple[i]
  pairDict = pairDictTuple[i]
  First = True
  for info in pair: # parse the specific INPUT_PAIRWISE_DE file
   if First == True: First = False # ignore the header
   else: 
    info = info.rstrip()
    info = info.split("\t")
    if gmt_file == True: pairDict[info[0]] = [info[4],info[7]] # for GO terms etc. from a gmt file
    else: pairDict[info[0]] = [info[6],info[10]] # for transcripts when using an xls file
 for GO in goDict: #iterate through relevant sequence IDs obtained from annotation file
  if gmt_file == True: print(GO + "\t" + goDict[GO], end="\t", file=OUTPUT_FILE) # for GO term annotations etc. from a gmt file
  # for transcript annotations when using an xls file from Trinotate report:
  else: print(GO + "\t" + goDict[GO][5] + "\t" + goDict[GO][8] + "\t" + goDict[GO][9] + "\t" + goDict[GO][10], end="\t", file=OUTPUT_FILE)
  for x in range(4): # iterate through each dictionary created for each INPUT_PAIRWISE_DE file
   pairDict = pairDictTuple[x]
   if GO in pairDict:
    print(pairDict[GO][0] + "\t" + pairDict[GO][1], end ="\t", file=OUTPUT_FILE)
   else: print("NA" + "\t" + "NA", end="\t", file=OUTPUT_FILE)
  print("\n", end="", file=OUTPUT_FILE)
