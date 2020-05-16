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

- `hapcounter.py` counts the total number of sequences and unique haplotypes per species and outputs to csv table
<details>
  <summary>hapcounter usage</summary>
  python hapcounter.py -h:

  usage: hapcounter.py [-h] [-i] [-f]

  optional arguments:

  -h, --help    show this help message and exit

  -i , --inputfile  Sequences input file name

  -f , --inputfileformat    'fasta' is recommended, other formats not tested
</details>

- `seq_purger.py` retains unique haplotypes only and removes sequences that do not meet minimum unambiguous length requirement

- Calculate all intraspecific distances in dataset (intra)
- Calculate all maximum intraspecific distances in dataset (Dmax)

- Calculate all interspecific distances in dataset (inter) [TO DO]
- Calculate all minimum distances to the nearest neighbor (Dmin_NN) [TO DO]

- the Jupyter Notebook `graphs.ipynb` contains scripts to interactively generate ('barcode gap') violin plots and output them for publication.

- `chao2.py` Uses all species' fastas in /species_fastas to run the SpideR_chao2.R script to calculate chao2 estimates of the total haplotype diversity and returns a csv. Note that chao2 assumes random sampling.

- `SpideR_haploaccum.R` R script that plots haplotype accumulation curves, based on the SpideR package (https://cran.r-project.org/web/packages/spider/spider.pdf)

# Example workflow

To do.
