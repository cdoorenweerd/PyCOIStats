from Bio import AlignIO
from collections import Counter
import argparse

# A site is parsimony-informative if it contains at least two types of nucleotides (or amino acids), 
# and at least two of them occur with a minimum frequency of two.

parser = argparse.ArgumentParser(description="Script to summarize number of variable and parsimony informative sites in an alignment")
parser.add_argument("-i", "--inputfile", metavar="", required=True, 
                    help="Sequences input file name")
parser.add_argument("-f", "--inputfileformat", metavar="", default='fasta',
                    help="'fasta' is recommended, other formats not tested")
args = parser.parse_args()

inputfile = args.inputfile
inputfileformat = args.inputfileformat

def aln_summarizer(inputfile, inputfileformat):
    sequences = AlignIO.read(inputfile, inputfileformat)
    totallength = sequences.get_alignment_length()
    print("records: " + str(len(sequences))) # this is the number of records not sequence length!
    print("alignment length: " + str(sequences.get_alignment_length())) # this is alignment length

    strict = ["A","C","T","G"]
    iupac = ["A","C","T","G","N","?","-","M","R","W","S","Y","K","V","H","D","B"]
    strictvariablesites = 0
    iupacvariablesites = 0
    strictparsimonysites = 0
    iupacparsimonysites = 0

    for pos in range(totallength):
        poslist = []
        for record in sequences:
            poslist.append(record.seq[pos]) # add all characters of position in all records to position list
        posdict = Counter(poslist) # dictionary with unique items of list as keys, values as number of times represented
        strictchars = [x for x in posdict.keys() if x in strict] # compare keylist against strict; all values must be present
        if len(posdict.keys()) > 1 and len(strictchars) > 1: # more than 1 key means its a variable site
            strictvariablesites += 1
            # need to count if there are two or more values > 1 to be parsimony informative
            if sum(value > 1 for value in posdict.values()) > 1:
                strictparsimonysites += 1
        iupacchars = [x for x in posdict.keys() if x in iupac]
        if len(posdict.keys()) > 1 and len(iupacchars) > 1:
            iupacvariablesites += 1
            if sum(value > 1 for value in posdict.values()) > 1:
                iupacparsimonysites += 1

    print("ACTG variable sites: " + str(strictvariablesites))
    print("ACTG parsimony informative sites: " + str(strictparsimonysites))
    print("IUPAC variable sites: " + str(iupacvariablesites))
    print("IUPAC parsimony informative sites: " + str(iupacparsimonysites))

aln_summarizer(inputfile,inputfileformat)
        
