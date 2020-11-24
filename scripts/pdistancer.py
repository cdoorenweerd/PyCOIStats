#/usr/bin/env python
# Camiel Doorenweerd 2020

import itertools
import pandas as pd
import argparse
import os
import time
from Bio import AlignIO
from basefunctions import IUPACdistance
from basefunctions import createlistofspecies


parser = argparse.ArgumentParser(description="Script to calculate all inter and intra pairwise distances, and intra Dmax and inter Dmin_NN and outputs to an overall csv and one with statistics per species.")
parser.add_argument("-i", "--inputfile", metavar="", required=True,
                    help="Alignment input file name")
parser.add_argument("-f", "--inputfileformat", metavar="", default='fasta',
                    help="'fasta' [default] is recommended, other formats not tested but may work")
args = parser.parse_args()


inputfile = args.inputfile
inputfileformat = args.inputfileformat
inputfileclean = os.path.splitext(inputfile)[0]
outputfile = str(str(inputfileclean) + "_pdistances.csv")
speciesstatsfile = str(str(inputfileclean) + "_speciesstats.csv")


def pcomparisons(inputfile, inputfileformat):
    starttime = time.perf_counter()
    sequences = AlignIO.read(inputfile, inputfileformat)
    print(sequences)
    print("Generating pairwise distance matrix, this may take a while...")
    pdist_dict = {}
    for a, b in itertools.combinations(sequences, 2):
        pdist = IUPACdistance(str(a.seq), str(b.seq))
        species1_id = str(a.id).split(".")[0]
        species1 = str(a.id).split(".")[1]
        species2_id = str(b.id).split(".")[0]
        species2 = str(b.id).split(".")[1]
        key = species1_id + "|" + species2_id
        pdist_dict.update({key: [species1_id, species1, species2_id, species2, pdist]})
    df_pdist = pd.DataFrame.from_dict(pdist_dict, orient='index', 
                                      columns=['species1_id', 'species1',
                                               'species2_id', 'species2',
                                               'pdist'])
    print(str(len(df_pdist.index)) + " pairwise comparisons")
    endtime = time.perf_counter()
    runtime = round((endtime - starttime),3)
    print(str(runtime) + " seconds to calculate pairwise distances")
    return df_pdist


def overallstats(df_pdist):
    starttime = time.perf_counter()
    print("Sorting intra- and interspecific values")
    df_intra = df_pdist[df_pdist['species1'] == df_pdist['species2']]
    print(str(len(df_intra.index)) + ' intraspecific distance values')
    df_inter = df_pdist[df_pdist['species1'] != df_pdist['species2']]
    print(str(len(df_inter.index)) + ' interspecific distance values')
    endtime = time.perf_counter()
    runtime = round((endtime - starttime),3)
    print(str(runtime) + " seconds to sort")
    return df_intra, df_inter


def neighbor_species(row,speciesname):
    if row['species1'] == speciesname:
        return row['species2']
    if row['species2'] == speciesname:
        return row['species1']


def speciesstats(df_pdist):
    starttime = time.perf_counter()
    print("Calculating statistics per species")
    listofspecies = createlistofspecies(inputfile, inputfileformat)
    dmaxvalues = []
    dmin_nnvalues = []
    sp_avg = {}
    for speciesname in listofspecies:
        df_intraperspecies = df_pdist[(df_pdist['species1'] == speciesname) & (df_pdist['species2'] == speciesname)]
        avg_intra = df_intraperspecies['pdist'].mean() # turn to N/A when empty?
        d_max = df_intraperspecies['pdist'].max() # turn to N/A when empty?
        dmaxvalues.append(d_max)
        n_intra_comparisons = len(df_intraperspecies.index)

        df_neighbors12 = df_pdist[(df_pdist['species1'] == speciesname) & (df_pdist['species2'] != speciesname)]
        df_neighbors21 = df_pdist[(df_pdist['species2'] == speciesname) & (df_pdist['species1'] != speciesname)]
        df_neighbors = pd.concat([df_neighbors12, df_neighbors21], ignore_index=True)
        df_neighbors['neighbor'] = df_neighbors.apply(lambda row: neighbor_species(row,speciesname), axis=1)
        avg_inter = df_neighbors['pdist'].mean()
        dmin_nn = df_neighbors['pdist'].min()
        dmin_nnvalues.append(dmin_nn)
        n_inter_comparisons = len(df_neighbors.index)
        nearestneighbor = df_neighbors[(df_neighbors['pdist'] == dmin_nn)]['neighbor'].iloc[0]

        sp_avg.update({speciesname: [avg_intra, d_max, n_intra_comparisons,
                                     avg_inter, dmin_nn, n_inter_comparisons,
                                     nearestneighbor]})
    df_sp_avg = pd.DataFrame.from_dict(sp_avg, orient='index', columns=['avg_intra', 'intra_d_max',
                                                                        'n_intra_comparisons',
                                                                        'avg_inter', 'inter_dmin_nn',
                                                                        'n_inter_comparisons',
                                                                        'nearest_neighbor'])
    df_intradmax = pd.DataFrame({'intra_dmax': dmaxvalues})
    df_interdmin_nn = pd.DataFrame({'inter_dmin_nn': dmin_nnvalues})
    endtime = time.perf_counter()
    runtime = round((endtime - starttime),3)
    print(str(runtime) + " seconds to calculate species statistics")
    return df_sp_avg, df_intradmax, df_interdmin_nn


df_pdist = pcomparisons(inputfile, inputfileformat)

df_intra, df_inter = overallstats(df_pdist)
df_sp_avg, df_intradmax, df_interdmin_nn = speciesstats(df_pdist)

df_sp_avg.to_csv(speciesstatsfile)
print("Species statistics written to " + str(speciesstatsfile))

df_distances = pd.concat([df_intra[['pdist']],df_intradmax,df_inter[['pdist']],df_interdmin_nn], ignore_index=False, axis=1) # need to clean up header names. # contains NULL values.
df_distances.to_csv(outputfile)
print("Overall intra, inter, Dmax and Dmin_NN distances written to " + str(outputfile))