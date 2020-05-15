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
parser.add_argument("-l", "--minlength", type=int, metavar="", default=0, 
                    help="Minimum length (bp) of sequences to be retained. Does not count N, ? or - towards length. Default '0'.")
args = parser.parse_args()


inputfile = args.inputfile
inputfileformat = args.inputfileformat
min_len = args.minlength
outputfile = str("qualifying_" + str(inputfile))


def seq_len_purge(inputfile, inputfileformat, min_len):
    min_len_seqs = []
    sequences = SeqIO.parse(inputfile, inputfileformat, alphabet=IUPAC.ambiguous_dna)
    for record in sequences:
        seq = str(record.seq)
        if (len(seq.translate(str.maketrans('','','N?-')))) > min_len:
            min_len_seqs.append(record)
    print("Sequences with min " + str(min_len) + " bp: " + str(len(min_len_seqs)))
    return min_len_seqs


def uniques(min_len_seqs, inputfile, inputfileformat):
    recordlist = []
    uniquedict = {}   
    for record in min_len_seqs:
        uniquedict[record.id] = record.seq
        recordlist.append(record)
    for a, b in itertools.combinations(recordlist, 2):
        seq1 = str(a.seq)
        seq2 = str(b.seq)
        if IUPACdistance(seq1, seq2) == 0:
            sp_seq1 = (a.id.split(".")[1])
            sp_seq2 = (b.id.split(".")[1])
            seq1_hq = (len(seq1.translate(str.maketrans('','','N?-MRWSYKVHDB'))))
            seq2_hq = (len(seq2.translate(str.maketrans('','','N?-MRWSYKVHDB'))))
            if b.id in uniquedict.keys() and sp_seq1 != sp_seq2:
                print("Warning: " + str(a.id) + " and " + str(b.id) + " share the same haplotype, both are retained.")
            elif b.id in uniquedict.keys():
                if seq1_hq > seq2_hq:
                    del uniquedict[b.id]
                else:
                    del uniquedict[a.id]
    print("Unique sequences (p-dist != 0) with minimum length: " + str(len(uniquedict))) 
    unique_seqs = []
    sequences = SeqIO.parse(inputfile, inputfileformat, alphabet=IUPAC.ambiguous_dna)
    for record in sequences:
        if record.id in uniquedict.keys():
            unique_seqs.append(record)
    return unique_seqs


count = 0
sequences = SeqIO.parse(inputfile, inputfileformat, alphabet=IUPAC.ambiguous_dna)
for record in sequences:
        count += 1
print("Found " + str(count) + " sequences in inputfile.")


min_len_seqs = seq_len_purge(inputfile, inputfileformat, min_len)
unique_min_len_seqs = uniques(min_len_seqs, inputfile, inputfileformat)


SeqIO.write(unique_min_len_seqs, outputfile, "fasta")
print("Qualifying sequences written to " + str(outputfile))