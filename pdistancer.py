#/usr/bin/env python
# Camiel Doorenweerd 2020

import itertools
import pandas as pd
import argparse
import os
from statistics import mean
from Bio import SeqIO
from basefunctions import IUPACdistance
from basefunctions import createlistofspecies


parser = argparse.ArgumentParser(description="Script to calculate all inter and intra pairwise distances, and intra Dmax and inter Dmin_NN and outputs to csv.")
parser.add_argument("-i", "--inputfile", metavar="", required=True,
                    help="Sequences input file name")
parser.add_argument("-f", "--inputfileformat", metavar="", default='fasta',
                    help="'fasta' [default] is recommended, other formats not tested")
args = parser.parse_args()


inputfile = args.inputfile
inputfileformat = args.inputfileformat
inputfileclean = os.path.splitext(inputfile)[0]
outputfile = str(str(inputfileclean) + "_pdistances.csv")
speciesstatsfile = str(str(inputfileclean) + "_speciesstats.csv")


def average(list):
    if len(list) > 0:
        avg = sum(list) / len(list)
    else:
        avg = 'N/A'
    return avg


pdistdict = []
sequences = SeqIO.parse(inputfile, inputfileformat)
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
        
listofspecies = createlistofspecies(inputfile, inputfileformat)
dmaxvalues = []
dmin_nnvalues = []
sp_avg = {}
for speciesname in listofspecies:
    intraperspecies = []
    interperspecies = []
    neighbors = {}
    for pair in pdistdict:
        species1 = str(pair).split(".")[1]
        species2 = str(pair).split(".")[3]
        pdist = float((str(pair).split(": ")[1]).replace("}", "")) # there is probably a better way to do this
        if speciesname == species1 == species2:
            intraperspecies.append(pdist)
        elif speciesname == species1 != species2:
            interperspecies.append(pdist)
            neighbors.update({species2: [pdist]})
        elif speciesname == species2 != species1:
            interperspecies.append(pdist)
            neighbors.update({species1: [pdist]})
    d_max = 'N/A'
    dmin_nn = 'N/A'
    nearestneighbor = 'N/A'
    d_nearestneighbor = 'N/A'
    sp_avg.update({speciesname: [average(intraperspecies),
                                 d_max,
                                 len(intraperspecies),
                                 average(interperspecies),
                                 dmin_nn,
                                 len(interperspecies),
                                 nearestneighbor,
                                 d_nearestneighbor]})
    if len(intraperspecies) > 0:
        d_max = max(intraperspecies)
        dmaxvalues.append(d_max)
        sp_avg[speciesname][1] = d_max
    if len(interperspecies) > 0:
        dmin_nn = min(interperspecies)
        dmin_nnvalues.append(dmin_nn)
        sp_avg[speciesname][4] = dmin_nn
    if len(neighbors) > 0:
        d_nearestneighbor = min(neighbors.values())
        sp_avg[speciesname][7] = d_nearestneighbor
        nearestneighbor = list(neighbors.keys())[list(neighbors.values()).index(d_nearestneighbor)]
        sp_avg[speciesname][6] = nearestneighbor



print(str(len(dmaxvalues)) + " intraspecific Dmax values.")
print(str(len(dmin_nnvalues)) + " interspecific Dmix_NN values.")

df_sp_avg = pd.DataFrame.from_dict(sp_avg, orient='index', columns=['avg_intra',
                                                                    'intra_d_max',
                                                                    'n_intra',
                                                                    'avg_inter',
                                                                    'inter_dmin_nn',
                                                                    'n_inter',
                                                                    'nearest_neighbor',
                                                                    'd_nearestneighbor'])
df_sp_avg.to_csv(speciesstatsfile)
print("P-distance averages per species written to " + str(outputfile))

df_intra = pd.DataFrame({'all_intra': intravalues})
df_intradmax = pd.DataFrame({'intra_dmax': dmaxvalues})
df_inter = pd.DataFrame({'all_inter': intervalues})
df_interdmin_nn = pd.DataFrame({'inter_dmin_nn': dmin_nnvalues})

df_distances = pd.concat([df_intra,df_intradmax,df_inter,df_interdmin_nn], ignore_index=False, axis=1)
df_distances.to_csv(outputfile)
print("All p-distances written to " + str(outputfile))