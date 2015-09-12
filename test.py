# Test the accuracy of the model. 
# Data was downloaded from PDBTM database, transmembrane domain and soluable domain was 
# extracted using parseTM.py and stored in pdbTM.dat and pdbNoneTM.dat
# Sensitivity is defined by the accurate positive prediction over the total count of 
# transmembrane domain. Specifity id defined by correct prediction over the predicted segments.
# The current model has sensitivity of 93.4%, and specifity of 88.5%
# Jiwei Xu@090815

from feature_generate import *


TMcorrect = 0
TMtotal = 0
with open('pdbTM.dat') as f:
    for line in f:
        line = line.strip()
        if len(line) >= 20:
            TMtotal += 1
            line = line[len(line)/2-10:len(line)/2+10]
            if predictFromSeq(line) > 0.5:
                TMcorrect += 1

noneTMtotal = 0
noneTMcorrect = 0
with open('pdbNoneTM.dat') as f:
    for line in f:
        line = line.strip()
        if line != '':
            noneTMtotal += 1
            if predictFromSeq(line) < 0.5:
                noneTMcorrect += 1


print 'Total peptides number: %d' % (TMtotal+noneTMtotal)
print 'Sensitivity is: %.3f' % (float(TMcorrect)/TMtotal)
print 'Specifity is: %.3f' % (float(TMcorrect)/(TMcorrect + noneTMtotal - noneTMcorrect))

