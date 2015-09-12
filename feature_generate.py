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

    :param seq (string): amino acids sequence in lenth of 20

    :return: nparray: of all features
             [Hydropathicity, side chain charge, polar, Interface Scale, pKa_sidechain, pI,
             Octanol Scale, number of hydrophobic resides, number of 'RK's]
    """
    features = [0]*19
    for aa in seq:
        for i in range(7):
            features[i] += aa_info[aa][i]
        for i in [0, 1, 2, 3]:
            if aa_info[aa][i] > 0:
                features[7+2*i] += aa_info[aa][i]
            else: 
                features[8+2*i] += aa_info[aa][i]
        if aa_info[aa][6] > 0:
            features[15] += aa_info[aa][i]
        else: 
            features[16] += aa_info[aa][i]
        if aa in 'VILMFWC':
            features[17] += 1
        if aa in 'RK':
            features[18] += 1
    return np.array(features)


def getP(x, w, intercept):
    """ calculate the probability of a set of features to be transmembrane domain

    :param x (nparray): all the features, nparray(9)
    :param w (nparray): model weight, nparray(9)
    :param intercept (float): model intercept

    :return: float: the probability of a set of features to be transmembrane domain
    """
    rawPrediction = x.dot(w)+intercept
    # Bound the raw prediction value
    rawPrediction = min(rawPrediction, 20)
    rawPrediction = max(rawPrediction, -20)
    return 1/(1+exp(-rawPrediction))


def predictFromSeq(seq):
    """ calculate the probability of a sequence to be transmembrane domain

    :param seq (string): amino acids sequences

    :return: float: the probability of a sequence to be transmembrane domain
    """
    # the value for w and intercept are from current best model
    w = np.array([0.127752353647,0.0248992252734,0.00743572396795,-0.103688386021,
                  -0.00750701102293,-0.00174062161425,-0.102516810939,-0.02997609561,
                  0.157728449257,-0.0161521904671,0.0410514157405,-0.0123707848711,
                  0.019806508839,-0.0484790847163,-0.0552093013047,-0.0478423792262,
                  -0.0558460067948,0.0278284300459,-0.0161521904671])
    intercept = 0.999916812177
    return getP(featureGenerate(seq), w, intercept)

if __name__ == '__main__':
    # test case
    seq0 = 'ERWAVLGGLLCLFMVRITMA'
    seq1 = 'WLMVKIFLILIVFFIIFVVFL'
    seq2 = 'MIPRKRYGSKNTDQGVYLGLSKTQVLSPATAGSSSSDIAPLPPPVTLVPPPPDTMSCRDRTQEFLSACKSLQTRQNGIQTNKPALRAVRQRSEFTLMAKRIGKDLSNTFAKLEKLTILAKRKSLFDDKAVEIEELTYIIKQDINSLNKQIAQLQDFVRAKGSQSGRHLQTHSNTIVVSLQSKLASMSNDFKSVLEVRTENLKQQRSRREQFSRAPVSALPLAPNHLGGGAVVLGAESHASKDVAIDMMDSRTSQQLQLIDEQDSYIQSRADTMQNIESTIVELGSIFQQLAHMVKEQEETIQRIDENVLGAQLDVEAAHSEILKYFQSVTSNRWLMVKIFLILIVFFIIFVVFLA'
    print predictFromSeq(seq0)
    print predictFromSeq(seq1[:20])
    for i in range(len(seq2)/20-1):
        print predictFromSeq(seq2[20*i:20*i+21])
