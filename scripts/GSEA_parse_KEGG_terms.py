#!/usr/bin/python3

'''
Title: GSEA_parse_KEGG_terms.py
Date: 2018-01-11
Author: Adam Nunn

Description:
    This program takes two input file: 1) a current "ko00000.keg" file, as
    downloaded from http://www.genome.jp/ containing KEGG terms and
    relevant metadata such as descriptions (INPUT_KEGG_TERMS), and 2) a
    previously constructed *.gmt format file for KEGG orthology terms and
    corresponding transcripts annotated with the term in the dataset
    (INPUT_GMT_FILE). The program will search the downloaded file for kegg
    pathway descriptions and corresponding kegg orthology terms, then
    combine all the transcripts annotated with relevant kegg orthology
    terms to generate a new *.gmt format file for kegg pathways containing
    pathway metadata and the corresponding transcripts (OUTPUT_GMT_FILE).

Usage:
    ./GSEA_parse_KEGG_terms.py INPUT_KEGG_TERMS INPUT_GMT_FILE OUTPUT_GMT_FILE
eg. ./GSEA_parse_KEGG_terms.py ko00000.keg KO.term.gmt KEGG.temp.gmt

'''

import sys

pathDict = {} # dictionary of each pathway and all corresponding orthology terms within them
descDict = {} # dictionary of each pathway and the corresponding description
koDict = {} # dictionary of each orthology term and all corresponding transcripts annotated with the term
flag = False # determines if pathway is relevant

with open(sys.argv[1], 'r') as INPUT_KEGG_TERMS, open(sys.argv[2], 'r') as INPUT_GMT_FILE, open(sys.argv[3], 'w') as OUTPUT_GMT_FILE:
 for info in INPUT_KEGG_TERMS: # parse the downloaded kegg terms to generate pathDict and descDict
  info = info.strip()
  info = info.split(" ")
  if (info[0] is "C") and info[4].startswith("0"): # find relevant kegg pathways
   if len(pathDict) > 1: pathDict[pathway] = set(pathDict[pathway]) # remove duplicates
   pathway = "ko" + info[4] # obtain kegg pathway ID
   desc = " ".join(info[5:-1]) # obtain corresponding pathway name/description
   pathDict[pathway] = [] # create entry for kegg pathway in pathDict, create empty list for orthology terms
   descDict[pathway] = desc # create entry for kegg pathway in descDict, give corresponding description
   flag = True # relevant pathway found!
  elif (info[0] is "C") and not info[4].startswith("0"): # filter out irrelevant pathway data
   if len(pathDict) > 1: pathDict[pathway] = set(pathDict[pathway]) # remove duplicates
   flag = False # pathway is irrelevant, prevent further orthology terms from appending to pathDict
  elif (info[0] is "D") and flag == True: # relevant pathway found, now append corresponding orthology terms to pathDict
   pathDict[pathway].append(info[6]) # add orthology term to the pathway entry in pathDict
 
 for line in INPUT_GMT_FILE: # parse orthology term GMT file to generate koDict
  line = line.strip()
  line = line.split("\t")
  koDict[line[0]] = line[2:]
 
 for x in enumerate(pathDict): # iterate through pathway terms in pathDict
  # print first two columns of GMT file output, check iteration number when formatting output
  if x[0] == 0: print(x[1] + "\t" + descDict[x[1]], end = "\t", file=OUTPUT_GMT_FILE)
  else: print("\n" + x[1] + "\t" + descDict[x[1]], end = "\t", file=OUTPUT_GMT_FILE)
  for i in pathDict[x[1]]: # iterate through each corresponding orthology term in pathDict[x[1]]
   if i in koDict: # check if specific orthology term exists in the dataset
    transcripts = koDict[i] # get all transcripts annotated with the term
    transcripts = "\t".join(transcripts) # join transcripts into a tab-delimited string
    print(transcripts, end = "\t", file=OUTPUT_GMT_FILE) # append all transcripts to output GMT file
