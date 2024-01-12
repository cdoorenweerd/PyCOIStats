#!/usr/bin/env python

import pandas as pd
import csv
import argparse
import os
import sys
from Bio import Phylo
sys.setrecursionlimit(10000) # without this script fails on large trees

parser = argparse.ArgumentParser(description="Test for monophyly of taxa in list in a Newick tree" )
parser.add_argument("-i", "--inputfile", metavar="", required=True,
                    help="Tree file in Newick format")
parser.add_argument("-t", "--taxalist", metavar="", required=True, 
                    help="Single column CSV list with taxa names to test for monophyly")
args = parser.parse_args()

inputtree = args.inputfile
taxa = args.taxalist # do we need to replace spaces with underscores?
inputtreeclean = os.path.splitext(inputtree)[0]
outputfile = str(str(inputtreeclean) + "_monophyly.csv")

def checkifintree(inputtree, flattaxalist):
    tree = Phylo.read(inputtree, "newick")
    Phylo.draw_ascii(tree)
    for taxon in flattaxalist:
        try:
            taxonstring = taxon + ".*" #regex matches species name and anything that follows
            tree.root_with_outgroup({"name": taxonstring}) # just to test if taxon is in the tree, will throw exception if not
        except Exception as e:
            print(taxon + " is not in this tree, skipping")
            print(e)
            flattaxalist.remove(taxon)
    return flattaxalist

def checkmonophyly(inputtree, flattaxalist):
    tree = Phylo.read(inputtree, "newick")
    tree.rooted = True
    speciesstatus = []
    for taxon in flattaxalist:
        taxonstring = taxon + ".*" #regex matches species name and anything that follows
        if tree.is_monophyletic(tree.find_clades({"name": taxonstring})) is not False:
            print(str(taxon) + " is monophyletic")
            speciesstatus.append({'taxon': taxon, 'monophyly': "monophyletic"})
        else:
            print(str(taxon) + " is non-monophyletic")
            speciesstatus.append({'taxon': taxon, 'monophyly': "non-monophyletic"})
    return speciesstatus

with open(taxa, "r", encoding="utf8") as taxafile:
    reader = csv.reader(taxafile)
    taxalist = list(reader)
    flattaxalist = [name for sublist in taxalist for name in sublist]

checkifintree(inputtree, flattaxalist)
if len(flattaxalist) == 0:
    print("None of the taxa were found in the tree; exiting")
    exit()

speciesstatus = checkmonophyly(inputtree, flattaxalist)
totaldf = pd.DataFrame(speciesstatus)
totaldf.to_csv(outputfile, index=False)
print("Output written to " + str(outputfile))