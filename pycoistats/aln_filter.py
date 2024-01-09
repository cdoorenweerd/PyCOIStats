#!/usr/bin/env python

from Bio import SeqIO
from Bio import AlignIO
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
parser.add_argument("-w", "--wobblestart", type=int, metavar="", default=0,
                    help="if not 0, removes the 3rd codon (wobble) positions from the alignment, using the number indicated here as the first codon position (1-3)")
args = parser.parse_args()

inputfile = args.inputfile
inputfileformat = args.inputfileformat
min_len = args.minlength
outputfile = str("q_" + str(inputfile))
wobblestart = args.wobblestart

def seq_len_filter(inputfile, inputfileformat, min_len):
    print("Removing sequences that do not meet the minimum length requirement")
    min_len_seqs = []
    sequences = AlignIO.read(inputfile, inputfileformat)
    for record in sequences:
        seq = str(record.seq)
        if (len(seq.translate(str.maketrans('','','N?-')))) > min_len:
            min_len_seqs.append(record)
    print("Sequences with min " + str(min_len) + " bp: " + str(len(min_len_seqs)))
    return min_len_seqs # this is a AlignIO object, like a dict with id, seq


def isdistinct(min_len_seqs, inputfile, inputfileformat):
    print("Removing non-distinct haplotypes across all samples by checking pairwise distance")
    recordlist = []
    specieslist = []
    distinctdict = {}
    for record in min_len_seqs:
        distinctdict[record.id] = record.seq # fills distinctdict with all min_len_seqs dict records
        recordlist.append(record)
    print("Starting with " + str(len(distinctdict)) + " sequences")
    totalcombinations = (len(list(itertools.permutations(recordlist, 2))))/2
    print("Calculating " + str(totalcombinations) + " pairwise distances")   
    progresscounter = 1
    for a, b in itertools.combinations(recordlist, 2):
        if a.id in distinctdict.keys() and b.id in distinctdict.keys():
            seq_a = str(a.seq)
            seq_b = str(b.seq)
            if IUPACdistance(seq_a, seq_b) == 0:
                sp_seq_a = a.id.split(".")[1]
                sp_seq_b = b.id.split(".")[1]
                if sp_seq_a != sp_seq_b:
                    print("Warning: " + str(a.id) + " and " + str(b.id) + " share the same haplotype, both are retained.")
                else:
                    seq_a_hq = len(seq_a.translate(str.maketrans('','','N?-MRWSYKVHDB')))
                    seq_b_hq = len(seq_a.translate(str.maketrans('','','N?-MRWSYKVHDB')))
                    if seq_a_hq > seq_b_hq:
                        del distinctdict[b.id]
                    else:
                        del distinctdict[a.id]
        progresscounter += 1
        percentcompletion = (progresscounter / totalcombinations) * 100
        print("Progress: " + str(round(percentcompletion, 3)) + " %                   ", end='\r')
    print(str(len(distinctdict)) + " sequences remain")
    print("Removing non-distinct haplotypes within species that were retained because they are shared across species")
    for record in min_len_seqs:
        species = (record.id.split(".")[1])
        specieslist.append(species)
    specieslist = list(set(specieslist)) # set() takes unique values only
    print("Found " + str(len(specieslist)) + " species")
    for species in specieslist:
        samples = []
        for record in min_len_seqs:
            if record.id.endswith(species):
                samples.append(record)
        for a, b in itertools.combinations(samples, 2):
            if a.id in distinctdict.keys() and b.id in distinctdict.keys():
                seq_a = str(a.seq)
                seq_b = str(b.seq)
                if IUPACdistance(seq_a, seq_b) == 0:
                    seq_a_hq = len(seq_a.translate(str.maketrans('','','N?-MRWSYKVHDB')))
                    seq_b_hq = len(seq_a.translate(str.maketrans('','','N?-MRWSYKVHDB')))
                    if seq_a_hq > seq_b_hq:
                        del distinctdict[b.id]
                    else:
                        del distinctdict[a.id]
    print("Final count of distinct sequences: " + str(len(distinctdict))) 
    distinct_seqs = []
    sequences = SeqIO.parse(inputfile, inputfileformat)
    for record in sequences:
        if record.id in distinctdict.keys():
            distinct_seqs.append(record)
    return distinct_seqs

def wobbleremover(wobblestart):
    startpos = (wobblestart-1) # python starts counting at 0
    parsed_file = "thirdcodonremoved_" + outputfile
    with open(outputfile) as original, open(parsed_file, 'w') as newfile:
        records = AlignIO.read(outputfile, fasta)
        for record in records:
            # retain first two characters of every set of three starting at startpos
            parsedrecord = record[startpos+1::3] 
            SeqIO.write(parsedrecord, newfile, 'fasta')


sequences = AlignIO.read(inputfile, inputfileformat)
print("Found " + str(len(sequences)) + " sequences in inputfile.")

min_len_seqs = seq_len_filter(inputfile, inputfileformat, min_len)
distinct_min_len_seqs = isdistinct(min_len_seqs, inputfile, inputfileformat)

SeqIO.write(distinct_min_len_seqs, outputfile, "fasta")
print("Qualifying sequences alignment written to " + str(outputfile))

if wobblestart > 0:
    wobbleremover(wobblestart)
    print("Qualifying sequences alignment without 3rd codon positions written to thirdcodonremoved_" + str(outputfile))
