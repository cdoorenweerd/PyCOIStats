#/usr/bin/env python
# Camiel Doorenweerd 2020


from Bio import SeqIO
from Bio.Alphabet import IUPAC
import itertools
import argparse
from basefunctions import IUPACdistance

parser = argparse.ArgumentParser(description="Script to remove duplicate haplotypes and select on minimum sequence bp length. Only works for protein coding sequences (gaps are regarded as missing data)" )
parser.add_argument("-i", "--inputfile", metavar="", 
                    help="Sequences input file name")
parser.add_argument("-f", "--inputfileformat", metavar="", default='fasta',
                    help="'fasta' is recommended, other formats not tested")
parser.add_argument("-l", "--minlength", metavar="", default=0, 
                    help="Minimum length (bp) of sequences to be retained. Does not count N, ? or - towards length. Default '0'.")
args = parser.parse_args()


inputfile = args.inputfile
inputfileformat = args.inputfileformat
min_len = args.minlength
outputfile = str("purged_" + str(inputfile))


def uniques(inputfile, inputfileformat):
    recordlist = []
    uniquedict = {}
    sequences = SeqIO.parse(inputfile, inputfileformat, alphabet=IUPAC.ambiguous_dna)
    for record in sequences:
        uniquedict[record.id] = record.seq
        recordlist.append(record)
    for a, b in itertools.combinations(recordlist, 2):
        if IUPACdistance(str(a.seq), str(b.seq)) == 0:
            sp_seq1 = (a.id.split(".")[1])
            sp_seq2 = (b.id.split(".")[1])
            if sp_seq1 != sp_seq2:
                print("Warning: " + str(a.id) + " and " + str(b.id) + " share the same haplotype.")
            if b.id in uniquedict.keys():
                del uniquedict[b.id]          
            
    print("Unique sequences (p-dist != 0): " + str(len(uniquedict)))
    
    unique_seqs = []
    sequences = SeqIO.parse(inputfile, inputfileformat, alphabet=IUPAC.ambiguous_dna)
    for record in sequences:
        if record.id in uniquedict.keys():
            unique_seqs.append(record)
    return unique_seqs


def seq_len_purge(unique_seqs, inputfileformat, min_len):
    unique_min_len_seqs = []
    for record in unique_seqs:
        seq = str(record.seq)
        if (len(seq.translate(str.maketrans('','','N?-')))) > min_len:
            unique_min_len_seqs.append(record)
    print("Unique sequences with min " + str(min_len) + " bp: " + str(len(unique_min_len_seqs)))
    return unique_min_len_seqs

count = 0
sequences = SeqIO.parse(inputfile, inputfileformat, alphabet=IUPAC.ambiguous_dna)
for record in sequences:
        count += 1
print("Found " + str(count) + " sequences in inputfile.")


unique_seqs = uniques(inputfile, inputfileformat)
unique_min_len_seqs = seq_len_purge(unique_seqs, inputfileformat, min_len)

SeqIO.write(unique_min_len_seqs, outputfile, "fasta")
print("Qualifying sequences written to " + str(outputfile))