#/usr/bin/env python
# Camiel Doorenweerd 2020

import itertools
import pandas as pd
import argparse
from basefunctions import IUPACdistance
from basefunctions import createlistofspecies


parser = argparse.ArgumentParser(description="Script to calculate all inter and intra pairwise distances, and intra Dmax and inter Dmin_NN and outputs to csv.")
parser.add_argument("-i", "--inputfile", metavar="", required=True,
                    help="Sequences input file name")
parser.add_argument("-f", "--inputfileformat", metavar="", default='fasta',
                    help="'fasta' is recommended, other formats not tested")
args = parser.parse_args()


inputfile = args.inputfile
inputfileformat = args.inputfileformat
outputfile = str(str(inputfile) + "_pdistances.csv")


pdistdict = []
sequences = SeqIO.parse(inputfile, inputfileformat, alphabet=IUPAC.ambiguous_dna)
for a, b in itertools.combinations(sequences, 2):
    pdist = IUPACdistance(str(a.seq), str(b.seq))
    pdistdict.append({(str(a.id) + '.' + str(b.id) + '.'): pdist})
print(str(len(pdistdict)) + " pairwise comparisons.")

intravalues = []
intervalues = []
for pair in pdistdict:
    species1 = str(pair).split(".")[1]
    species2 = str(pair).split(".")[3]
    pdist = float((str(pair).split(": ")[1]).replace("}", "")) # there is probably a better way to do this
    if species1 == species2:
        intravalues.append(pdist)
    else:
        intervalues.append(pdist)
print(str(len(intravalues)) + " intraspecific values.")
print(str(len(intervalues)) + " interspecific values.")
        
listofspecies = createlistofspecies(inputfile, fileformat)
dmaxvalues = []
dmin_nnvalues = []
for speciesname in listofspecies:
    intraperspecies = []
    interperspecies = []
    for pair in pdistdict:
        species1 = str(pair).split(".")[1]
        species2 = str(pair).split(".")[3]
        pdist = float((str(pair).split(": ")[1]).replace("}", "")) # there is probably a better way to do this
        if speciesname == species1 == species2:
            intraperspecies.append(pdist)
        elif speciesname == species1 != species2:
            interperspecies.append(pdist)
    if len(intraperspecies) > 0:
        dmaxvalues.append(max(intraperspecies))
    if len(interperspecies) > 0:
        dmin_nnvalues.append(min(interperspecies)) 
print(str(len(dmaxvalues)) + " intraspecific Dmax values.")
print(str(len(dmin_nnvalues)) + " interspecific Dmix_NN values.")

df_intra = pd.DataFrame({'all_intra': intravalues})
df_intradmax = pd.DataFrame({'intra_dmax': dmaxvalues})
df_inter = pd.DataFrame({'all_inter': intervalues})
df_interdmin_nn = pd.DataFrame({'inter_dmin_nn': dmin_nnvalues})

df_distances = pd.concat([df_intra,df_intradmax,df_inter,df_interdmin_nn], ignore_index=False, axis=1)
df_distances.to_csv(outputfile)
print("Results written to " + str(outputfile))