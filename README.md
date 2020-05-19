[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![Anaconda](https://img.shields.io/conda/v/fastai/fastai.svg)](https://anaconda.org/fastai/fastai)

# PyCOIStats package
A collection of scripts to manage COI data, calculate pairwise distances and plot them

### Dependencies

All dependencies are available through Anaconda (https://anaconda.org/anaconda/repo). It is recommended to create a new conda environment with the required packages:

- pandas
- itertools
- biopython
- numpy
- argparse
- os
- matplotlib
- seaborn
- csvkit
- r-base

e.g. using

```conda create -n pycoistats python=3.6 pandas itertools biopython numpy argpase os matplotlib seaborn csvkit r-base```

See further instructions on how to use conda environments at https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html


### Unique haplotypes

A unique haplotype is defined as a confidently different sequence; i.e. IUPAC ambiguity codes or missing data ("?") and gaps ("-") are ignored in comparisons.

For example:

AACTGTCA and AACTNY-A

Are considered 'identical' haplotypes (or you might say 'compatible' or 'non-unique').


The scripts compare all sequences in a pairwise fashion, so the computational time increases exponentially with more sequences. However, it should be able to handle ~5,000 sequences in less than 30 minutes on most desktop machines, for larger datasets a cluster is advisable.


### Fasta input format

ALL SEQUENCES MUST BE PROPERLY ALIGNED

The scripts were created for COI sequence data but will work on any protein coding aligned mitochondrial DNA sequence dataset. The functions will work on nuclear DNA as well, but ambiguity codes there can mean a) uncertainty in the data or b) different alleles -- it makes more sense to count ambiguities as differences with nuDNA.

### Species name recognition

The scripts will use anything after a period "." in the sample identifier as the species name. Sample identifiers must be unique for each sequence. So for example this would work:

>'>LEASV1523-19.Caloptilia_betulicola
NACTCTTTATTTTATTTTTGGAATTTGAGCCGGTATATTAGGAACTTCTTTAAGAATATTAATTCGAGCAGAATTAGGTAATCCAGGATCTTTAATTGGGGATGATCAAATTTATAATACAATTGTTACAGCTCATGCTTTCATTATAATTTTCTTTATAGTTATACCTATTATAATTGGGGGATTTGGGAATTGATTAGTCCCATTGATATTAGGAGCACCTGATATAGCTTTCCC
>
>'>ABASV1522-19.Caloptilia_rufipennella
NACTCTTTATTTTATTTTTGGAATTTGATCCGGTATATTAGGAACTTCTTTAAGAATATTAATTCGAGCAGAGTTAGGTAATCCAGGATCTTTAATTGGTGATGATCAAATTTATAATACCATTGTTACAGCTCACGCTTTTATTATAATTTTTTTTATAGTTATACCTATTATAATTGGGGGATTTGGAAATTGATTAGTGCCATTAATATTAGGGGCACCTGATATAGCATTCCC
>
>'>LETRA472-19.Caloptilia_rufipennella
NACTCTCTACTTTATTTTCGGAATTTGATCTGGAATATTAGGAACATCTTTAAGTATATTAATTCGAGCTGAATTAGGTAATCCAGGATCTTTAATTGGGGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTTTATAGTTATACCTATTATAATTGGTGGATTTGGAAACTGATTAGTGCCATTAATATTAGGGGCTCCTGATATAGCTTTCCC


### Script functions

For each script, run python script.py -h for usage instructions.


- `hapcounter.py` counts the total number of sequences and unique haplotypes per species and outputs to csv table.

- `seq_purger.py` retains sequences that meet the minimum length requirement and are unique haplotypes, and writes the result to a new fasta. Identical haplotypes across different species issue a warning and are retained.

- `pdistancer.py` calculates all intraspecific distances (all_intra), all interspecific distances (all_inter), the maximum intraspecific distances (Dmax) and minimum distances to the nearest neighbor (Dmin_NN) and outputs to a csv.

- the Jupyter Notebook `graphs.ipynb` contains scripts to interactively generate ('barcode gap') violin plots from the output from ```pdistancer.py``` and output the graphs for publication.

- `chao1.py` Uses all species' fastas in /species_fastas to run the SpideR_chao1.R script to calculate chao 1 estimates of the total haplotype diversity and returns a csv. Note that the function assumes a large number of specimens have been sampled and that duplicate haplotypes have not been removed.

- `SpideR_haploaccum.R` R script that plots haplotype accumulation curves, based on the SpideR package (https://cran.r-project.org/web/packages/spider/spider.pdf)

# Example workflow

To do.
