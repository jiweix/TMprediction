
# coding: utf-8

# This is try to predict transembrane domaim from large set of data
# 

# In[2]:

from pyspark.mllib.regression import LabeledPoint
from pyspark.mllib.classification import LogisticRegressionWithSGD
from math import log, exp


# In[3]:

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


# In[4]:

def parse(line):
    """ Extract labels and features from raw data
    
    :param line: single line from the input file, starts with 0 or 1, 
                 0 means that it is not transmenbrane domain, 1 means it is transmenbrane domain
                 20 charactors after the space, sanding for 20 residues
        
    :return: LabeledPoint: labeled with 0 or 1, and the features calculated from the peptides sequence  
             [Hydropathicity, side chain charge, polar, MW, pKa_sidechain, pI]
             
    """
    allAminoAcids = 'ACDEFGHIKLMNPQRSTVWY'
    label, seq = line.split()
    features = [0]*9
    for aa in seq:
        for i in range(7):
            features[i] += aa_info[aa][i]
        if aa in 'VILMFWC':
            features[7] += 1
        if aa in 'RK':
            features[8] += 1
    return LabeledPoint(label, features)


# In[5]:

fileName = 'YeastTM20.dat'
rawData = sc.textFile(fileName, 2).map(parse)
print rawData.take(5)


# In[6]:

weights = [0.8, 0.1, 0.1]
seed = 1
rawTrainData, rawValidationData, rawTestData = rawData.randomSplit(weights, seed)
rawTrainData.cache()
rawValidationData.cache()
rawTestData.cache()
nTrain = rawTrainData.count()
nVal = rawValidationData.count()
nTest = rawTestData.count()
print nTrain, nVal, nTest, nTrain + nVal + nTest


# Loss should calculated for a give prediction and label

# In[7]:

def computeLogLoss(p, y):
    """Calculates the value of log loss for a given probabilty and label.

    Note:
        log(0) is undefined, so when p is 0 we need to add a small value (epsilon) to it
        and when p is 1 we need to subtract a small value (epsilon) from it.

    :param p (float): A probabilty between 0 and 1.
           y (int): A label.  Takes on the values 0 and 1.

    :return: float: The log loss value.
    """
    epsilon = 10e-12
    if y == 1:
        pp = p
    if y == 0:
        pp = 1-p
    if pp == 0:
        return -log(pp+epsilon)
    elif pp ==1:
        return -log(pp-epsilon)
    else:
        return -log(pp)
    
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


# In[8]:

def evaluateResults(model, data):
    """ Calculates the log loss for the data given the model.

    :param model (LogisticRegressionModel): A trained logistic regression model.
           data (RDD of LabeledPoint): Labels and features for each observation.

    :return: float: Log loss for the data.
    """
    log_loss = (data.map(lambda x: computeLogLoss(getP(x.features, model.weights, model.intercept), x.label))
                    .reduce(lambda x, y: x+y))/data.count()
    return log_loss


# In[10]:

# try fixed hyperparameters
numIters = 500
stepSize = 1
regParam = 1e-6
regType = 'l2'
includeIntercept = True

model0 = LogisticRegressionWithSGD.train(rawTrainData,
                                         iterations=numIters, 
                                         step=stepSize, 
                                         miniBatchFraction=1.0, 
                                         initialWeights=None, 
                                         regParam=regParam, 
                                         regType=regType, 
                                         intercept=includeIntercept)
print model0.weights, model0.intercept


# In[11]:

classOneFracTrain = (rawTrainData.map(lambda x: x.label)
                                 .reduce(lambda x, y: x+y))/rawTrainData.count()
print classOneFracTrain

logLossTrBase = (rawTrainData.map(lambda x: x.label)
                             .map(lambda x: computeLogLoss(classOneFracTrain, x))
                             .reduce(lambda x, y: x+y))/rawTrainData.count()
print 'Baseline Train Logloss = {0:.3f}\n'.format(logLossTrBase)


# In[12]:

logLossTrLR0 = evaluateResults(model0, rawTrainData)
print ('Logloss:\n\tLogReg = {0:.3f}'
       .format(logLossTrLR0))


# In[13]:

logLossVa = evaluateResults(model0, rawValidationData)
print ('Logloss:\n\tLogReg = {0:.3f}'
       .format(logLossVa))


# In[15]:

numIters = 100
regType = 'l2'
includeIntercept = True

# Initialize variables using values from initial model training
bestModel = None
bestLogLoss = 1e10

stepSizes = [0.01, 0.1, 1, 10]
regParams = [1e-6, 1e-3]
for stepSize in stepSizes:
    for regParam in regParams:
        model = (LogisticRegressionWithSGD
                 .train(rawTrainData, numIters, stepSize, regParam=regParam, regType=regType,
                        intercept=includeIntercept))
        logLossVa = evaluateResults(model, rawValidationData)
        print ('\tstepSize = {0:.2f}, regParam = {1:.0e}: logloss = {2:.3f}'
               .format(stepSize, regParam, logLossVa))
        if (logLossVa < bestLogLoss):
            bestModel = model
            bestLogLoss = logLossVa

print ('Validation Logloss:\n\tBaseline = {0:.3f}\n\tLogReg = {1:.3f}'
       .format(logLossTrBase, bestLogLoss))


# In[16]:

logLossTe = evaluateResults(bestModel, rawTestData)
print ('Logloss:\n\tLogReg = {0:.3f}'
       .format(logLossTe))


# In[10]:

# More iteration with optimized parameters
numIters = 5000
stepSize = 0.01
regParam = 1e-6
regType = 'l2'
includeIntercept = True

model1 = LogisticRegressionWithSGD.train(rawTrainData,
                                         iterations=numIters, 
                                         step=stepSize, 
                                         miniBatchFraction=1.0, 
                                         initialWeights=None, 
                                         regParam=regParam, 
                                         regType=regType, 
                                         intercept=includeIntercept)
print model1.weights, model1.intercept


# In[18]:

logLossTe1 = evaluateResults(model1, rawTestData)
print ('Logloss:\n\tLogReg = {0:.3f}'
       .format(logLossTe1))


# In[28]:

for x in rawData.take(10):
    print x
    print getP(x.features, model1.weights, model1.intercept)
    print computeLogLoss(getP(x.features, model.weights, model.intercept), x.label)


# In[29]:

print model0.weights


# In[ ]:



