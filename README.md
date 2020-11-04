[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![DOI](https://zenodo.org/badge/264048060.svg)](https://zenodo.org/badge/latestdoi/264048060)

# PyCOIStats package
A collection of scripts to analyze COI data, filter non-distinct haplotypes, calculate pairwise distances and plot graphs of this data.

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


### Distinct haplotypes

A principle difference between PyCOIStats and most other software is in how it defines a distinct haplotype. A distinct haplotype is a confidently different sequence; i.e. IUPAC ambiguity codes or missing data ("?") and gaps ("-") are ignored in comparisons.

For example for deciding whether two haplotypes are distinct:

AACTGTCA and AACTNY-A

are considered 'identical' haplotypes (or you might say 'compatible' or 'non-distinct').

This principle also affects the pairwise distance estimates. Assume for example the following two sequences:

GTAAYTNN and GTAACTGC

Most other software counts 1 difference in 8 bases: 12.5% difference. PyCOIStats counts 1 difference in 6 bases (ignoring the ambiguities and not counting them towards compared strand length): 16.7% difference

Ideally, when COI data is of extremely high quality and all sequences are of equal length and have no ambiguous bases, this would not make a difference. However, practice is that missing or ambiguous data is common and should be interpreted correctly.


### Computational demand

The scripts compare all sequences in a pairwise fashion, so the computational time increases exponentially with more sequences. However, it should be able to handle ~5,000 sequences in a couple of hours on most desktop machines, for larger datasets a computing cluster is advisable.


### Fasta input format

ALL SEQUENCES MUST BE PROPERLY ALIGNED

The scripts were created for COI sequence data but will work on any protein coding aligned mitochondrial DNA sequence dataset. The functions will work on nuclear DNA as well, but ambiguity codes there can mean a) uncertainty in the data or b) different alleles -- it makes more sense to count ambiguities as differences with nuDNA.

### Species name recognition

SAMPLE NAMING CONVENTION

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

For each script, run `python script.py -h` for usage instructions.


- `hapcounter.py` counts the total number of sequences and distinct haplotypes per species and outputs to csv file.

- `seq_purger.py` filters sequences to meet the minimum length requirement and be distinct haplotypes, and writes the result to a new fasta. Non-distinct haplotypes across different species raise a warning and are retained.

- `pdistancer.py` calculates all intraspecific distances (all_intra), all interspecific distances (all_inter), the maximum intraspecific distances (Dmax) and minimum distances to the nearest neighbor (Dmin_NN) and outputs to a csv. Per v1.2 and later, a second csv table is produced with statistics per species [intra_Dmax, n_intra comparisons, avg_inter distance, inter_Dmin_nn, n inter comparison, nearest neighbor].


- the Jupyter Notebook `graphs.ipynb` contains scripts to interactively generate ('barcode gap') violin plots from the csv output from ```pdistancer.py``` and output the graphs for publication.

- `makespeciesfastas.py` will generate a separate fasta for each species in the folder /species_fastas

- `chao1.py` Uses all species' fastas in /species_fastas to run the SpideR_chao1.R script to calculate chao 1 estimates of the total haplotype diversity and returns a csv. Note that the function assumes a large number of specimens have been sampled and that duplicate haplotypes have not been removed.

- `SpideR_haploaccum.R` R script that plots haplotype accumulation curves, based on the SpideR package (https://cran.r-project.org/web/packages/spider/spider.pdf)
