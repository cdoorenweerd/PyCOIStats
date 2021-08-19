#!/usr/bin/env python
# Camiel Doorenweerd 2020
# For a full list of alignment formats that it can handle see https://biopython.org/wiki/SeqIO

import argparse
import pandas as pd
import os
from Bio import AlignIO
from Bio import SeqIO


parser = argparse.ArgumentParser(description="Rename sequence identifiers using a reference CSV list")
parser.add_argument("-i", "--inputfile", metavar="", required=True,
                    help="FASTA input file")
parser.add_argument("-r", "--renamelist", metavar="", 
                    help="CSV file with first column with header 'current_name' and corresponding new names in a second column named 'new_name'.")
parser.add_argument("-f", "--format", metavar="", default ='fasta',
                    help="Input file format, such as fasta, nexus, phylip. [default = fasta]")
parser.add_argument("-l", "--listcurrentnames", action="store_true",
                    help="If this flag is given, script will only produce a 'current_names.txt' file with all current sequence names")
args = parser.parse_args()


inputfile = args.inputfile
inputfileformat = args.format
outputfile = os.path.splitext(inputfile)[0] + "_renamed.fas"
listrequest = args.listcurrentnames


if listrequest == True:
    with open('current_names.txt', 'w') as currentnameslist:
        aln = AlignIO.read(inputfile, inputfileformat)
        for record in aln:
            print(record.id, file=currentnameslist)
else:
    df = pd.read_csv(args.renamelist, header=0)
    df.set_index('current_name', inplace=True)
    with open(outputfile, 'w') as replaced:
        aln = AlignIO.read(inputfile, inputfileformat)
        print("Found " + str(len(aln)) + " sequences in input file")
        count = 0
        for record in aln:   
            try:
                newname = df.loc[record.id, 'new_name']
                record.id = newname
                count +=1
            except KeyError:
                record.id = record.id # don't change name if not in csv list
                print("Could not match " + record.id)
            SeqIO.write(record, replaced, 'fasta-2line')
        print("Successfully changed " + str(count) + " records and written to " + outputfile)