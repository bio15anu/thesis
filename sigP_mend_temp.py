#!/usr/bin/python3

'''
Title: sigP_mend_temp.py
Date: 2018-01-11
Author: Adam Nunn

Description:
    This program works alongside "sigP_mend_Transpep_for_sigP.py" to fix
    file formating issues when running SignalP on sequences with long
    sequence IDs. The program takes two input files: 1) the output file
    generated by SignalP after successful completion of analysis
    (INPUT_SIGP_FILE), and 2) the original fasta file that was to be
    analysed by SignalP in the first place (INPUT_FASTA_FILE) NOTE: this
    file must be unchanged since running "sigP_mend_Transpep_for_sigP.py".
    The program will reformat the file generated by SignalP so that it
    will correspond with the sequence IDs in INPUT_FASTA_FILE, generating
    a new SignalP output file (OUTPUT_SIGP_FILE).

Usage:
    ./sigP_mend_temp.py INPUT_SIGP_FILE INPUT_FASTA_FILE OUTPUT_SIGP_FILE
eg. ./sigP_mend_temp.py signalp.temp Trinity_isoH_cdhit10.fasta.transdecoder.pep signalp.out

'''

import sys

count = 0
idList = []
countList = []

with open(sys.argv[1], 'r') as INPUT_SIGP_FILE, open(sys.argv[2], 'r') as INPUT_FASTA_FILE, open(sys.argv[3], 'w') as OUTPUT_SIGP_FILE:
 for line in INPUT_FASTA_FILE: # parse through the protein fasta file containing original sequence IDs
  line = line.rstrip()
  if line.startswith(">"):
   line = line.split(" ")[0]
   idList.append(line[1:])
   countList.append("ID_{}".format(count))
   count += 1
  else: continue
 mydict = dict(zip(countList, idList))
 for info in INPUT_SIGP_FILE: # parse through the file generated by SignalP, replace sequence IDs with original IDs
  info = info.rstrip()
  if info.startswith("ID"):
   info = info.split("\t")
   ID = info[0]
   stringy = '\t'.join(info[1:])
   print(mydict[ID], end="\t", file=OUTPUT_SIGP_FILE)
   print(stringy, file=OUTPUT_SIGP_FILE)
  else: print(info, file=OUTPUT_SIGP_FILE)
   
