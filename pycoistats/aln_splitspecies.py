#!/usr/bin/env python
# Camiel Doorenweerd 2020

import argparse
from Bio import AlignIO
from pathlib import Path
from basefunctions import createlistofspecies

parser = argparse.ArgumentParser(description="Script to create a fasta per species in a subfolder.")
parser.add_argument("-i", "--inputfile", metavar="", 
                    help="Sequences input file name")
parser.add_argument("-f", "--inputfileformat", metavar="", default='fasta',
                    help="'fasta' [default] is recommended, other formats not tested")
parser.add_argument("-m", "--min_seqs", type=int, metavar="", default=0,
                    help="Minimum number of sequences to write a fasta for that species, default 0")
args = parser.parse_args()


inputfile = args.inputfile
inputfileformat = args.inputfileformat
min_seqs = args.min_seqs

Path("./species_fastas").mkdir(parents=True, exist_ok=True)

listofspecies = createlistofspecies(inputfile, inputfileformat)

for speciesname in listofspecies:
    recordlist = []
    sequences = AlignIO.read(inputfile, inputfileformat)
    for record in sequences:
        if (record.id.split(".")[1]) == speciesname:        
            recordlist.append(record)
    print(str(len(recordlist)) + " sequence(s) for " + str(speciesname))
    if len(recordlist) >= min_seqs:
        AlignIO.write(recordlist, ("./species_fastas/" + str(speciesname) + ".fas"), "fasta")