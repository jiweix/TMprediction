# TMprediction
# Update 082915

The predictFromSeq function in feature_generate.py could be used to predict whether a 20 residues peptide is transmenbrane domain. 
The return value of that function is the probability of a sequence to be transmembrane domain

Current model is calculated from 19 features generated from a sequence:
[Hydropathicity, side chain charge, polar, Interface Scale, pKa_sidechain, pI, Octanol Scale, 
Hydropathicity positive only, Hydropathicity negative only, side chain charge positive only, side chain charge negative only,
polar positive only, polar negative only, Interface Scale positive only, Interface Scale negative only, 
Octanol Scale positve only, Octanol Scale negative only, number of hydrophobic resides, number of 'RK's]

TMprediction.ipynb used spark logistic regression model in mllib to calculate the model, raw data was divided into train, 
validation and test data. Different parameters of logistic regression were tried, stepSize 0.01 and regParam 1e-6 were used to 
generate a final model
