#!/usr/bin/env Rscript

library("ape")
library("spider")
library("optparse")
 
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


aln <- read.dna(opt$file, format='fasta')
chaoHaplo(aln)