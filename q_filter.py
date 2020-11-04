#/usr/bin/env python
# Camiel Doorenweerd 2020


from Bio import SeqIO
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


def seq_len_filter(inputfile, inputfileformat, min_len):
    min_len_seqs = []
    sequences = SeqIO.parse(inputfile, inputfileformat)
    for record in sequences:
        seq = str(record.seq)
        if (len(seq.translate(str.maketrans('','','N?-')))) > min_len:
            min_len_seqs.append(record)
    print("Sequences with min " + str(min_len) + " bp: " + str(len(min_len_seqs)))
    return min_len_seqs


def isdistinct(min_len_seqs, inputfile, inputfileformat):
    recordlist = []
    distinctdict = {}   
    for record in min_len_seqs:
        distinctdict[record.id] = record.seq
        recordlist.append(record)
    for a, b in itertools.combinations(recordlist, 2):
        if a.id in distinctdict.keys() and b.id in distinctdict.keys():
            seq_a = str(a.seq)
            seq_b = str(b.seq)
            if IUPACdistance(seq_a, seq_b) == 0:
                sp_seq_a = (a.id.split(".")[1])
                sp_seq_b = (b.id.split(".")[1])
                if sp_seq_a != sp_seq_b:
                    print("Warning: " + str(a.id) + " and " + str(b.id) + " share the same haplotype, both are retained.")
                else:
                    seq_a_hq = (len(seq_a.translate(str.maketrans('','','N?-MRWSYKVHDB'))))
                    seq_b_hq = (len(seq_a.translate(str.maketrans('','','N?-MRWSYKVHDB'))))
                    if seq_a_hq > seq_b_hq:
                        del distinctdict[b.id]
                    else:
                        del distinctdict[a.id]
    print("Distinct sequences (p-dist != 0) with minimum length: " + str(len(distinctdict))) 
    distinct_seqs = []
    sequences = SeqIO.parse(inputfile, inputfileformat)
    for record in sequences:
        if record.id in distinctdict.keys():
            distinct_seqs.append(record)
    return distinct_seqs


count = 0
sequences = SeqIO.parse(inputfile, inputfileformat)
for record in sequences:
        count += 1
print("Found " + str(count) + " sequences in inputfile.")


min_len_seqs = seq_len_filter(inputfile, inputfileformat, min_len)
distinct_min_len_seqs = isdistinct(min_len_seqs, inputfile, inputfileformat)


SeqIO.write(distinct_min_len_seqs, outputfile, "fasta")
print("Qualifying sequences written to " + str(outputfile))