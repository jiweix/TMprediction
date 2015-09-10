# TMprediction
# Update 082915

As current transmembrane domain prediction tools are not quite accurate, I tried to solve this by applying machine learning methods. 
I download transmembrane domain data from http://wardlab.cbs.umn.edu/yeast/ for the generation of the model, all transmembrane helixes with length larger than 20 was added to the database. To further validate my model, I download and extract the transmembrane domain from PDBTM (http://pdbtm.enzim.hu/). This data is very accurate because it is based on experimental data from X-ray crystallography, literally without error. 

Different parameters were used for transmembrane domain prediction, and they all yielded good results. Thus, I think the combination of these parameters would likely give better result. Current model is calculated from 19 features generated from a sequence:
[Hydropathicity, side chain charge, polar, Interface Scale, pKa_sidechain, pI, Octanol Scale, 
Hydropathicity positive only, Hydropathicity negative only, side chain charge positive only, side chain charge negative only,
polar positive only, polar negative only, Interface Scale positive only, Interface Scale negative only, 
Octanol Scale positve only, Octanol Scale negative only, number of hydrophobic resides, number of 'RK's]
These features covers most of the important parameters of a peptides, and the model calculated from these features was very accurate. 

The predictFromSeq function in feature_generate.py could be used to predict whether a 20 residues peptide is transmenbrane domain. The return value of that function is the probability of a sequence to be transmembrane domain

The model was test against the experimental data from PDBTM, and it accurately predict 93.5% transmembrane domain in the lenth range of 20-26 residues, which means this model has very high sensitivity. According to Y. Diao et al, the sensitivity of HMMTOP was 74.9% and that of TMpred was 66.4%, which are much lower than the current model in this project. The specifity of this model, defining by correct prediction over the predicted segments, is 88.5%, which is also higher than the two tools mentioned above. 

Ref:
Using pseudo amino acid composition to predict transmembrane regions in protein: celluar automata and Lempel-Ziv complexity. Y.Diao et al, Amino Acids (2008) 34: 111-117
