#/usr/bin/env python
# Camiel Doorenweerd 2020

from Bio import SeqIO
from Bio.Alphabet import IUPAC
import itertools
import os
import argparse
import pandas as pd


parser = argparse.ArgumentParser(description="Script to count number of sequences and unique haplotypes per species. Species names are parsed based on '.'")
parser.add_argument("-i", "--inputfile", metavar="", required=True, 
                    help="Sequences input file name")
parser.add_argument("-f", "--inputfileformat", metavar="", default='fasta',
                    help="'fasta' is recommended, other formats not tested")
args = parser.parse_args()


inputfile = args.inputfile
inputfileformat = args.inputfileformat
inputfileclean = os.path.splitext(inputfile)[0]
outputfile = str(str(inputfileclean) + "_haplotypecounts.csv")


def IUPACdistance(seq1, seq2):
    ignorelist = ["N","?","-","M","R","W","S","Y","K","V","H","D","B"]
    unamblengthseq1 = len(seq1.translate(str.maketrans('','','N?-MRWSYKVHDB')))
    unamblengthseq2 = len(seq2.translate(str.maketrans('','','N?-MRWSYKVHDB')))
    comparedlength = min(unamblengthseq1, unamblengthseq2)
    difference = 0
    for x, y in zip(seq1.upper(), seq2.upper()):
        if x in ignorelist or y in ignorelist:
            difference += 0
        elif x != y:
            difference += 1
    dpairwise = (difference / comparedlength)
    return dpairwise


def createlistofspecies(inputfile, fileformat):
    sequences = SeqIO.parse(inputfile, fileformat, alphabet=IUPAC.ambiguous_dna)
    listofspecies = []
    for record in sequences:
        speciesname = (record.id.split(".")[1])
        if speciesname not in listofspecies:
            listofspecies.append(speciesname)
    print(str(len(listofspecies)) + " species in charactermatrix.")
    return listofspecies


def hapcounter(listofspecies, inputfile, inputfileformat):
    speciesstats = []
    for speciesname in listofspecies:
        recordlist = []
        sequences = SeqIO.parse(inputfile, inputfileformat, alphabet=IUPAC.ambiguous_dna)
        for record in sequences:
            if (record.id.split(".")[1]) == speciesname:        
                recordlist.append(record)
        uniquelist=[]
        for record in recordlist:
            uniquelist.append(record.id)
        for a, b in itertools.combinations(recordlist, 2):
            if IUPACdistance(str(a.seq), str(b.seq)) == 0:
                if b.id in uniquelist:
                    uniquelist.remove(b.id)
        speciesstats.append({'species': speciesname,
                             'n_seq': len(recordlist),
                             'n_haplotypes': len(uniquelist)})
    return speciesstats

listofspecies = createlistofspecies(inputfile, inputfileformat)
speciesstats = hapcounter(listofspecies, inputfile, inputfileformat)

totaldf = pd.DataFrame(speciesstats)
totaldf.to_csv(outputfile, index=False)
print("Output written to " + str(outputfile))