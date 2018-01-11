#!/usr/bin/python3

'''
Title: GSEA_parse_GO_terms.py
Date: 2018-01-11
Author: Adam Nunn

Description:
    This program takes two input file: 1) a current "go.obo" file, as
    downloaded from http://www.geneontology.org/ containing GO terms and
    relevant metadata such as descriptions (INPUT_GO_TERMS), and 2) a
    partially constructed *.gmt format file where the column containing
    GO term descriptions is currently missing (INPUT_GMT_FILE). The
    program will search the database for current descriptions and append
    the information to the GMT file in a new output (OUTPUT_GMT_FILE).

Usage:
    ./GSEA_parse_GO_terms.py INPUT_GO_TERMS INPUT_GMT_FILE OUTPUT_GMT_FILE
eg. ./GSEA_parse_GO_terms.py go.obo GO.temp.gmt GO.term.gmt

'''

import sys

goid = []
godesc = []
flag = False

with open(sys.argv[1], 'r') as INPUT_GO_TERMS, open(sys.argv[2], 'r') as INPUT_GMT_FILE, open(sys.argv[3], 'w') as OUTPUT_GMT_FILE:
 for info in INPUT_GO_TERMS: # parse the downloaded file containing GO term metadata
  info = info.strip()
  if info.startswith("id:") and (flag == False):
   goid.append(info[4:])
   flag = True
  elif info.startswith("name:") and (flag == True):
   name = info[6:]
   godesc.append(name)
   flag = False
  elif info.startswith("alt_id:") and (flag == False):
   goid.append(info[8:])
   godesc.append(name)
 godict = dict(zip(goid, godesc))
 for line in INPUT_GMT_FILE: # parse the partially constructed GMT file
  line = line.split("\t")
  line.insert(1, godict[line[0]])
  line = '\t'.join(line)
  print(line.rstrip(), file=OUTPUT_GMT_FILE) # create new GMT file with GO term descriptions
