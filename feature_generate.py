from math import exp, log
import numpy as np
# A dictionary for parameters of single amino acid
# 'A': [Hydropathicity, side chain charge, polar, Interface Scale, pKa_sidechain, pI, Octanol Scale]
aa_info = {'A': [1.80, 0, 0, 0.17, 0.0, 6.0, 0.5],
           'C': [2.50, 0, 0, -0.24, 8.3, 5.1, 0],
           'D': [-3.5, -1, 1, 1.23, 3.9, 2.8, 3.64],
           'E': [-3.5, -1, 1, 0.0, 4.3, 3.2, 0.11],
           'F': [2.80, 0, 0, -1.13, 0.0, 5.5, -1.71],
           'G': [-0.4, 0, 0, 0, 0.0, 6.0, 1.15],
           'H': [-3.2, 0, -1, 0.17, 6.0, 7.6, 0.11],
           'I': [4.50, 0, 0, -0.31, 0.0, 6.0, -1.12],
           'K': [-3.9, 1, -1, 0.99, 10.5, 9.7, 2.8],
           'L': [3.80, 0, 0, -0.56, 0.0, 6.0, -1.25],
           'M': [1.90, 0, 0, -0.23, 0.0, 5.7, -0.67],
           'N': [-3.5, 0, 1, 0.42, 0.0, 5.4, 0.85],
           'P': [2.80, 0, 0, 0.45, 0.0, 6.3, 0.14],
           'Q': [-3.5, 0, 1, 0.58, 0.0, 5.7, 0.77],
           'R': [-4.5, 1, -1, 0.81, 12.5, 10.8, 1.81],
           'S': [-0.8, 0, 1, 0.13, 0.0, 5.6, 0.46],
           'T': [-0.7, 0, 1, 0.14, 0.0, 5.6, 0.25],
           'V': [4.20, 0, 0, 0.07, 0.0, 6.0, -0.46],
           'W': [-0.9, 0, 0, -1.85, 0.0, 5.9, -2.09],
           'Y': [-1.3, 0, 1, -0.94, 10.7, 5.7, -0.71]}


def featureGenerate(seq):
    """ generate a nparray features from amino acids sequence

    :param seq: amino acids sequence in lenth of 20

    :return: nparray: of all features
             [Hydropathicity, side chain charge, polar, Interface Scale, pKa_sidechain, pI,
             Octanol Scale, number of hydrophobic resides, number of 'RK's]
    """
    features = [0]*9
    if len(seq) > 20:
        seq = seq[:20]
    if len(seq) < 20:
        return features
    for aa in seq:
        for i in range(7):
            features[i] += aa_info[aa][i]
        if aa in 'VILMFWC':
            features[7] += 1
        if aa in 'RK':
            features[8] += 1
    return np.array(features)


def getP(x, w, intercept):
    """ calculate the probability of a set of features to be transmembrane domain

    :param x: all the features, nparray(9)
    :param w: model weight, nparray(9)
    :param intercept: model intercept

    :return: float: the probability of a set of features to be transmembrane domain
    """
    rawPrediction = x.dot(w)+intercept
    # Bound the raw prediction value
    rawPrediction = min(rawPrediction, 20)
    rawPrediction = max(rawPrediction, -20)
    return 1/(1+exp(-rawPrediction))


def predictFromSeq(seq):
    """ calculate the probability of a sequence to be transmembrane domain

    :param seq: amino acids sequences

    :return: float: the probability of a sequence to be transmembrane domain
    """
    # the value for w and intercept are from current best model
    w = np.array([0.179153018013, 0.0283969783774, -0.000504497958365, -0.120192164047,
                  -0.0173080116018, -0.037859079514, -0.131153338053, 0.0228327769761, -0.0181199689548])
    intercept = 0.999521602093
    return getP(featureGenerate(seq), w, intercept)

if __name__ == '__main__':
    # test case
    seq1 = 'IVLCSSGTLICCENAVVVLII'
    print predictFromSeq(seq1[:20])