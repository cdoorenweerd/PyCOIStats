#/usr/bin/env python
# Camiel Doorenweerd 2020

from Bio import AlignIO
import itertools
import os
import argparse
import pandas as pd
from skbio.diversity import alpha_diversity
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
    sequences = AlignIO.read(inputfile, inputfileformat)
    speciesstats = []
    progresscounter = 1
    for speciesname in listofspecies:
        print('Calculating stats for species ' + str(progresscounter) + ': ' + speciesname + '               ', end='\r')
        recordlist = []
        for record in sequences:
            if (record.id.split(".")[1]) == speciesname:        
                recordlist.append(record)
        uniquedict={}
        for record in recordlist:
            uniquedict.update({record.id: 1})
        for a, b in itertools.combinations(recordlist, 2):
            if IUPACdistance(str(a.seq), str(b.seq)) == 0:
                newvalue = uniquedict[a.id] + uniquedict[b.id]
                uniquedict.update({a.id: newvalue})
                uniquedict.update({b.id: 0})
        uniquedict = {key:value for key, value in uniquedict.items() if value != 0}
        hapcounts = list(uniquedict.values())
        chao1 = alpha_diversity('chao1', hapcounts)
        chao1_ci = alpha_diversity('chao1_ci', hapcounts)
        speciesstats.append({'species': speciesname,
                             'n_seq': len(recordlist),
                             'n_haplotypes': len(uniquedict),
                             'chao1': chao1[0],
                             'chao1_95ci': chao1_ci[0]
                             })
        progresscounter += 1
    print('')
    return speciesstats

listofspecies = createlistofspecies(inputfile, inputfileformat)
speciesstats = hapcounter(listofspecies, inputfile, inputfileformat)

totaldf = pd.DataFrame(speciesstats)
totaldf.to_csv(outputfile, index=False)
print("Output written to " + str(outputfile))