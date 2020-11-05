#/usr/bin/env python
#Camiel Doorenweerd 2018

import argparse
import pandas as pd
from Bio import SeqIO


parser = argparse.ArgumentParser(description="A script to remove the 3rd codon positions from a FASTA file" )
parser.add_argument("-i", "--inputfile", metavar="", required=True,
                    help="FASTA input file")
parser.add_argument("-s", "--startpos", metavar="", required=True, type=int, 
                    help="Start of first codon position in the alignment")
args = parser.parse_args()

startpos = (args.startpos-1) # python starts counting at 0
original_file = args.inputfile
parsed_file = "thirdcodonremoved_" + original_file

with open(original_file) as original, open(parsed_file, 'w') as newfile:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        # retain first two characters of every set of three starting at startpos
        parsedrecord = record[startpos+1::3] 
        SeqIO.write(parsedrecord, newfile, 'fasta')