__author__ = 'KATRINA'

'''

python d3pipeline.py [input DF file] [input classification file] [config file]
same as the original d3pipeline_1.py, but this one transforms the dataframe into a binary frame prior to apply tsne...attempting to see if that will solve the variance issue

#outdated

'''


import pandas as pd
#import sklearn.feature_selection as fs
import sklearn
from sklearn import svm, cross_validation, feature_selection
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.pipeline import make_pipeline
from sklearn.pipeline import Pipeline
import matplotlib as plt
import numpy as np
from sklearn.decomposition import PCA
import os
import json
import math
from sklearn import (manifold, datasets, decomposition, ensemble, random_projection)
import sys


def removeAllZeros(df):
    df_no_zeros = df.loc[(df!=0).any(axis=1)] #remove all zero-only rows
    return df_no_zeros

def removeRows(df, threshold):
    df_no_under_ten = df.loc[(df.sum(axis=1) >= threshold)] #remove all rows that have sum < 10
    return df_no_under_ten

def removeCols(df,threshold):
    df_no_col_under_twothousand = df.loc[:,(df.sum(axis=0) >= threshold)]
    return df_no_col_under_twothousand


def normalizeRPM(df):
    total_transcripts_per_cell = []
    for i in df.columns.values:
        sum = 0
        for genes in df[i]:
            sum += genes
        total_transcripts_per_cell.append(sum)

    total_sum = np.sum(total_transcripts_per_cell)
    per_million_scaling_factor = total_sum/1000000
    new_df = df.apply(lambda x: (x/per_million_scaling_factor))
    return new_df


input_file = sys.argv[1]
classification_file = sys.argv[2]

df = pd.DataFrame.from_csv((input_file),sep="\t")
classification_vector = json.load(open(sys.argv[2]))["Diagnosis"]
print(classification_vector)

df_classification = pd.DataFrame({'classif': classification_vector})#.transpose()
#x = [df.transpose(), df_classification.transpose()]
df.loc[len(df)]=classification_vector
df_with_class = df #pd.concat(x)

print(df_with_class)
print(df_with_class.iloc[0:len(df_with_class)-2])
print(df_with_class.iloc[len(df_with_class)-1])


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
else:
    gene_set_method = 'kbest'
    print("Default gene set method used: kbest")

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
print(df_clean)
x = df_clean.iloc[0:len(df_clean)-1]#df_clean.index != 'classif']#df_clean.drop(df.index['classif'])
classification_vector = df_clean.iloc[len(df_clean)-1]#[df_clean.index == 'classif']
df_clean = x
print(df_clean)
print(len(df_clean.columns.values))
print(len(classification_vector))
#df_clean = df #TODO: THIS IS TEMPORARY, BUT I NEED A WAY TO REMOVE COLUMNS FROM DF AND THEN REMOVE CORRESPONDING CLASSIFICATION VALUES
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

print(genes)
subset_df = df_clean[df_clean.index.isin(genes)]
subset_df.to_csv(os.path.join(os.path.dirname(sys.argv[1]),"_DFresult_bin.txt"),sep="\t")

variance = subset_df.var(axis=0) #variance in columns

if dim_red_method == 'TSNE':
    X = subset_df.transpose()
    n_samples, n_features = X.shape[0],X.shape[1]
    print(n_samples)
    print(n_features)
    print("begin tSNE")
    tsne = manifold.TSNE(n_components=3, early_exaggeration=2.0, init='pca', random_state=0)
    X_tsne = tsne.fit_transform(X)
    # plot the result
    vis_x = X_tsne[:,0]
    vis_y = X_tsne[:, 2] #1]
    print("completed tSNE")
    print(len(vis_x))
    print(len(vis_y))
    print(len(variance))
    print(len(classification_vector.as_matrix()))
    df1 = pd.DataFrame({'tSNEx': vis_x, 'tSNEy':vis_y, 'variance':variance,'classif':classification_vector.as_matrix()})

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
(output_df.transpose()).to_csv(os.path.join(os.path.dirname(sys.argv[1]),"d3data_bin.tsv"),sep="\t") #THIS FILE CAN BE LOADED INTO TSNE ON D3.js


#TODO: THIS DOES NOT INCLUDE THE NORAMLIZATION STEP!!
