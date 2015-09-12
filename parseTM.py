# This script extract transmembrane domain and soluable domain from data downloaded from PDBTM
# Data was stored in file 'pdbtmalpha.xml'
# Jiwei Xu@090815

TM = []
notTM = []

def valid(seq):
    """ test whether a given sequence is amino acids sequence

    :param seq (string): a given string to test

    :return: boolean: True if the string is a valid amino acids sequence
    """
    for char in seq:
        if char not in 'ACDEFGHIKLMNPQRSTVWY':
            return False
    return True


def extractTM(seq):
    """ find all transmembrane domain in the sequence, according to the notation in the xml file, then
        add these transmembrane domain into TM.
        If no TM in sequence, then divide the sequence into peptides with length of 20, and add these 
        peptides into notTM.

    :param seq (string): string extract from the xml file containing the sequence and 
                    TM information
    """
    seqStart = seq.find('<SEQ>')+5
    seqEnd = seq.find('</SEQ>')
    proteinSeq = seq[seqStart:seqEnd].strip().replace(' ', '').replace('\n', '')
    if not valid(proteinSeq):
        return
    find = False
    for line in seq[seqEnd+6:].split('\n'):
        if 'type="H"' in line:
            start = int(line[line.find('seq_beg="')+9:line.find('" pdb_beg')])
            end = int(line[line.find('seq_end="')+9:line.find('" pdb_end')])
            if end-start > 19 and end-start < 27:
                if proteinSeq[start:end+1] not in TM:
                    TM.append(proteinSeq[start:end+1])
            find = True
    if not find:
        for i in range(len(proteinSeq)/20):
            notTM.append(proteinSeq[i*20:(i+1)*20])


# read from file, get the information with label CHAIN in the xml file, and process the information with 
# extractTM
with open('pdbtmalpha.xml') as f:
    seq = ''
    toAdd = False
    for line in f:
        if '<CHAIN' in line:
            toAdd = True
        if '</CHAIN>' in line:
            toAdd = False
            extractTM(seq)
            seq = ''
        if toAdd:
            seq += line

# write output file
with open('pdbTM.dat', 'w') as f:
    for line in TM:
        f.write('%s\n' % line.encode('utf8'))

with open('pdbNoneTM.dat', 'w') as f:
    for line in notTM:
        f.write('%s\n' % line.encode('utf8'))  
