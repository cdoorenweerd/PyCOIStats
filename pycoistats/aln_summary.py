#!/usr/bin/env python

from Bio import AlignIO
from collections import Counter
import argparse

# A site is parsimony-informative if it contains at least two types of nucleotides (or amino acids), 
# and at least two of them occur with a minimum frequency of two.

parser = argparse.ArgumentParser(description="Script to summarize number of variable and parsimony informative sites in an alignment")
parser.add_argument("-i", "--inputfile", metavar="", required=True, 
                    help="Sequences input file name")
parser.add_argument("-f", "--inputfileformat", metavar="", default='fasta',
                    help="'fasta' is recommended, other formats not tested")
args = parser.parse_args()

inputfile = args.inputfile
inputfileformat = args.inputfileformat

def aln_summarizer(inputfile, inputfileformat):
    sequences = AlignIO.read(inputfile, inputfileformat)
    totallength = sequences.get_alignment_length()
    print(f"{'records: ' : <40}{len(sequences) :^5}")
    print(f"{'alignment length: ' : <40}{sequences.get_alignment_length() : ^5}")

    missing = ["N","?","-"]
    strict = ["A","C","T","G"]
    iupac = ["A","C","T","G","M","R","W","S","Y","K","V","H","D","B"]
    strictvariablesites = 0
    iupacvariablesites = 0
    strictparsimonysites = 0
    iupacparsimonysites = 0
    missingdatasites = 0

    for pos in range(totallength):
        poslist = []
        for record in sequences:
            poslist.append(record.seq[pos]) # add all characters of position in all records to position list
        posdict = Counter(poslist) # dictionary with unique items of list as keys, values as number of times represented
        # ACTG count
        strictchars = [x for x in posdict.keys() if x in strict] # compare keylist against strict; all values must be present
        if len(posdict.keys()) > 1 and len(strictchars) > 1: # more than 1 key means its a variable site
            strictvariablesites += 1
            # need to count if there are two or more values > 1 to be parsimony informative
            if sum(value > 1 for value in posdict.values()) > 1:
                strictparsimonysites += 1
        # IUPAC count
        iupacchars = [x for x in posdict.keys() if x in iupac]
        if len(posdict.keys()) > 1 and len(iupacchars) > 1:
            iupacvariablesites += 1
            if sum(value > 1 for value in posdict.values()) > 1:
                iupacparsimonysites += 1
        # missing count
        if posdict['N'] + posdict['?'] + posdict['-'] == len(sequences):
            missingdatasites += 1

    print(f"{'Sites with only missing data: ' : <40}{missingdatasites : ^5}{round(missingdatasites / totallength, 4) : ^15}")
    print(f"{'ACTG variable sites: ' : <40}{strictvariablesites : ^5}{round(strictvariablesites / totallength, 4) : ^15}")
    print(f"{'ACTG parsimony informative sites: ' : <40}{strictparsimonysites : ^5}{round(strictparsimonysites / totallength, 4) : ^15}")
    print(f"{'IUPAC variable sites: ' : <40}{iupacvariablesites : ^5}{round(iupacvariablesites / totallength, 4) : ^15}")
    print(f"{'IUPAC parsimony informative sites: ' : <40}{iupacparsimonysites : ^5}{round(iupacparsimonysites / totallength, 4) : ^15}")

def missingoverall(inputfile, inputfileformat):
    sequences = AlignIO.read(inputfile, inputfileformat)
    totalbp = sequences.get_alignment_length() * len(sequences)
    missingcount = 0
    for record in sequences:
        missing = record.seq.count("N") + record.seq.count("?") + record.seq.count("-")
        missingcount += missing
    missingoverall = missingcount / totalbp
    print(f"{'Missing overall: ' : <40}{missingcount : ^5}{round(missingoverall, 4) : ^15}")

print(f"{'' : <40}{'count' : ^5}{'proportion' : ^15}")
aln_summarizer(inputfile, inputfileformat)
missingoverall(inputfile, inputfileformat)
