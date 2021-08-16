from Bio import AlignIO


def IUPACdistance(seq1, seq2):
    ignorelist = ["N","?","-","M","R","W","S","Y","K","V","H","D","B"]
    unamblengthseq1 = len(seq1.translate(str.maketrans('','','N?-MRWSYKVHDB')))
    unamblengthseq2 = len(seq2.translate(str.maketrans('','','N?-MRWSYKVHDB')))
    comparedlength = min(unamblengthseq1, unamblengthseq2)
    difference = 0
    for x, y in zip(seq1.upper(), seq2.upper()):
        if x in ignorelist or y in ignorelist:
            difference += 0
        elif x != y:
            difference += 1
    dpairwise = (difference / comparedlength)
    return dpairwise


def createlistofspecies(inputfile, fileformat):
    sequences = AlignIO.read(inputfile, fileformat)
    listofspecies = []
    for record in sequences:
        speciesname = (record.id.split(".")[1])
        if speciesname not in listofspecies:
            listofspecies.append(speciesname)
    print(str(len(listofspecies)) + " species in charactermatrix.")
    return listofspecies