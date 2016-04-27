__author__ = 'KATRINA'

'''

from a counts file, classification file, and config file (specifying the method to use
for gene selection, thresholds, and dimensionality reduction), create a d3data.tsv file which can then be
used with the DropSeqViewer. Also creates an _DFresult.txt file containing only the subset
of genes selected.

python d3pipeline.py [input DF file] [input classification file] [config file]

example config file:
{
"gene_set_method":"kbest",
"gene_number":20,
"cell_expression_threshold":2000,
"dim_red_method":"ZIFA"
}

note:
- input dataframe is expected to have cells in the column orientation (and genes in rows)
- limited error checking on the config file, if you are missing a field in the config file
the script will likely error. TODO
'''


import pandas as pd
from sklearn import feature_selection
from ZIFA import ZIFA
import numpy as np
from sklearn.decomposition import PCA
import os
import json
from sklearn import (manifold)
import sys

def removeCols(df,threshold):
    df_no_col_under_twothousand = df.loc[:,(df.sum(axis=0) >= threshold)]
    return df_no_col_under_twothousand

input_file = sys.argv[1]
classification_file = sys.argv[2]

df = pd.DataFrame.from_csv((input_file),sep="\t")
classification_vector = json.load(open(sys.argv[2]))["classification"]
print(classification_vector)

df_classification = pd.DataFrame({'classif': classification_vector})#.transpose()
#x = [df.transpose(), df_classification.transpose()]
df.loc[len(df)]=classification_vector
df_with_class = df

config_file = json.load(open(sys.argv[3]))

dim_red_method = 'TSNE'
try:
    dim_red_method = config_file["dim_red_method"]
except:
    print("Default dimensionality reduction method used: TSNE")

genes = []
gene_number = 500

if config_file["gene_set_method"]:
    gene_set_method = config_file["gene_set_method"]

if gene_set_method == 'kbest':
    if config_file["gene_set_method"]:
        gene_number = config_file["gene_number"]
    else:
        gene_number = 500   #TODO: SET THIS DEFAULT VALUE BASED ON THE RESULTS OF THE EXPERIMENT FROM LupusTest
        print("Default gene number being used: 500")
elif gene_set_method == 'manual':
    if config_file["genes"]:
        genes = config_file["genes"]
    else:
        genes = []
        print("WARNING: Manual gene set was selected, but genes are not specified; will run kbest to select genes")

if "gene_expression_threshold" in config_file.keys():
    gene_expression_threshold = int(config_file["gene_expression_threshold"])
else:
    gene_expression_threshold = 10

if "cell_expression_threshold" in config_file.keys():
    cell_expression_threshold = int(config_file["cell_expression_threshold"])
else:
    cell_expression_threshold = 3000


#df_clean = removeCols(removeRows(df,gene_expression_threshold),cell_expression_threshold)
df_clean = removeCols(df_with_class,cell_expression_threshold)
x = df_clean.iloc[0:len(df_clean)-1]#df_clean.index != 'classif']#df_clean.drop(df.index['classif'])
classification_vector = df_clean.iloc[len(df_clean)-1]#[df_clean.index == 'classif']
df_clean = x
df_trans = df_clean.transpose()

if len(genes) == 0: # select k best
    X_new = feature_selection.SelectKBest(feature_selection.chi2, k=gene_number).fit_transform(df_trans, classification_vector)
    n = feature_selection.SelectKBest(feature_selection.chi2,k=gene_number).fit(df_trans, classification_vector)

    t = n.get_support()
    indices_of_interest = [i for i, x in enumerate(t) if x]
    genes_of_interest = []

    for i in indices_of_interest:
        genes_of_interest.append(df_trans.columns.values[i])

    genes = genes_of_interest

subset_df = df_clean[df_clean.index.isin(genes)]
subset_df.to_csv(os.path.join(os.path.dirname(sys.argv[1]),"_DFresult.txt"),sep="\t")

variance = subset_df.var(axis=0) #variance in columns

if dim_red_method == 'ZIFA':
    f = lambda x: np.log(1+x)
    logDF = subset_df.applymap(f) # DF_final.applymap(f)

    transposed_ZIFA = logDF.transpose()
    Z_trans, MP_trans = ZIFA.fitModel(transposed_ZIFA.as_matrix(),2)

    X=[]
    Y=[]
    for i in Z_trans:
        X.append(i[0])
        Y.append(i[1])
    df1 = pd.DataFrame({'tSNEx': X, 'tSNEy':Y, 'variance':variance,'classif':classification_vector.as_matrix()})

if dim_red_method == 'TSNE':
    pca = PCA(n_components=15)
    pcaF = pca.fit(df_clean) #THIS WILL NOT TAKE A SUBSET OF THE GENES

    X= (pca.components_).transpose()
    #X = subset_df.transpose()
    n_samples, n_features = X.shape[0],X.shape[1]
    print(n_samples)
    print(n_features)
    print("begin tSNE")
    tsne = manifold.TSNE(n_components=2, init='pca', random_state=0)#,perplexity=15)
    X_tsne = tsne.fit_transform(X)
    # plot the result
    vis_x = X_tsne[:,0]
    vis_y = X_tsne[:, 1] #1]
    print("completed tSNE")
    df1 = pd.DataFrame({'tSNEx': vis_x, 'tSNEy':vis_y, 'variance':variance,'classif':classification_vector.as_matrix()})
    subset_df = df_clean #TEMPORARY FIX TO MAKE IT SO THAT TSNE DOESN"T REMOVE ALL GENES FROM THE FILE

if dim_red_method == 'PCA':
    print("begin PCA")
    X = subset_df.transpose()
    pca = PCA(n_components=2)
    pca_res = pca.fit_transform(X)
    pca_x = pca_res[:,0]
    pca_y = pca_res[:,1]
    print("completed PCA")
    df1 = pd.DataFrame({'tSNEx': pca_x, 'tSNEy':pca_y, 'variance':variance,'classif':classification_vector.as_matrix()})

#df1 = pd.DataFrame({'tSNEx': vis_x, 'tSNEy':vis_y, 'variance':variance,'classif':classification_vector})
output_df = pd.concat([subset_df,df1.transpose()])
(output_df.transpose()).to_csv(os.path.join(os.path.dirname(sys.argv[1]),"d3data.tsv"),sep="\t") #THIS FILE CAN BE LOADED INTO TSNE ON D3.js




