#/usr/bin/env python
# Camiel Doorenweerd 2020

from Bio import AlignIO
import itertools
import os
import argparse
import pandas as pd
from basefunctions import IUPACdistance
from basefunctions import createlistofspecies


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


def hapcounter(listofspecies, inputfile, inputfileformat):
    speciesstats = []
    for speciesname in listofspecies:
        recordlist = []
        sequences = AlignIO.read(inputfile, inputfileformat)
        for record in sequences:
            if (record.id.split(".")[1]) == speciesname:        
                recordlist.append(record)
        uniquelist=[] # change to dict to allow for counts for chao1 estimates
        for record in recordlist:
            uniquelist.append(record.id) # add all records as keys and value as 1
        for a, b in itertools.combinations(recordlist, 2):
            if IUPACdistance(str(a.seq), str(b.seq)) == 0:
                if b.id in uniquelist:
                    uniquelist.remove(b.id)
                    # uniquedict(a.id) value += 1
        # create haplotype count list as 'counts = [3, 1, 4]'
        # calculate chao 1: 
        # from skbio.diversity import alpha_diversity
        # skbio.diversity.alpha.chao1(counts)
        # skbio.diversity.alpha.chao1_ci(counts) gives 95% lower bound, upper bound
        speciesstats.append({'species': speciesname,
                             'n_seq': len(recordlist),
                             'n_haplotypes': len(uniquelist)})
    return speciesstats

listofspecies = createlistofspecies(inputfile, inputfileformat)
speciesstats = hapcounter(listofspecies, inputfile, inputfileformat)

totaldf = pd.DataFrame(speciesstats)
totaldf.to_csv(outputfile, index=False)
print("Output written to " + str(outputfile))