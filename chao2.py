#/usr/bin/env python
# Camiel Doorenweerd 2020

import os
import subprocess
import pandas as pd


listoffastas = os.listdir('./species_fastas/')
print("Calculating Chao2 for " + str(len(listoffastas)) + " species")
outputfile = 'chao2_estimates.csv'


speciesstats = []
for fasta in listoffastas:
    speciesname = os.path.splitext(fasta)[0]
    print("Running " + str(speciesname))
    cmd = str("Rscript --vanilla SpideR_chao2.R -f ./species_fastas/") + fasta
    chao2result = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    speciesstats.append({'species': speciesname,
                         'chao2result':chao2result})

totaldf = pd.DataFrame(speciesstats)
totaldf.to_csv(outputfile, index=False)
print("Output written to " + str(outputfile))