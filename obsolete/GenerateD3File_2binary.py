__author__ = 'KATRINA'

'''

python GenerateD3File_2.py [input dataframe] [input classification .json file] [-o=output_file_name] [-s=[GENE1,GENE2,GENE3]]
same as the original GenerateD3File_1.py, but this one transforms the dataframe into a binary frame prior to apply tsne...attempting to see if that will solve the variance issue

note: this was an early version preceeding d3pipeline_2bin.py

#outdated

'''

import pandas as pd
import os
import json
from sklearn import (manifold, datasets, decomposition, ensemble, random_projection)
import sys

print("begin")
#GET DATA AND CLASSIFICATION VECTOR SOMEHOW:
#prefix = "C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\d3_test"
prefix = "C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\Lupus_data\\analysis_1"
input_file = sys.argv[1]
#df = pd.DataFrame.from_csv(os.path.join(prefix,"1000_df_result.tsv"),sep="\t")
df = pd.DataFrame.from_csv((input_file),sep="\t")

classification_vector = json.load(open(sys.argv[2]))["Diagnosis"]

if '-s' in sys.argv:
    genes_of_interest = sys.argv[4].split('=')[1].replace('[','').replace(']','').split(',')
    subset_df = df[df.index.isin(genes_of_interest)]
    subset_df.to_csv(os.path.dirname(sys.argv[1])+sys.argv[3].split('=')[1] + "_df_bin.tsv",sep="\t")
    df = subset_df

print("loaded data and classification")

N = df.transpose()
print(N)
N[N != 0] = 1 #convert to binary matrix
X = N
#X = df.transpose()

print(X)
n_samples, n_features = X.shape[0],X.shape[1]
print(n_samples)
print(n_features)
print("begin tSNE")
tsne = manifold.TSNE(n_components=2, early_exaggeration=32.0, init='pca', random_state=0)
X_tsne = tsne.fit_transform(X)
# plot the result
vis_x = X_tsne[:, 0]
vis_y = X_tsne[:, 1]
print("completed tSNE")


variance = df.var(axis=0) #variance in columns

#df1 = pd.DataFrame({'tSNEx': vis_x, 'tSNEy':vis_y, 'variance':variance,'classif':classification_vector},columns = [df.columns.values])

print(len(vis_x))
print(len(variance))
print(len(classification_vector))
df1 = pd.DataFrame({'tSNEx': vis_x, 'tSNEy':vis_y, 'variance':variance,'classif':classification_vector})#,index=['tSNEx','tSNEy','variance', 'classif'])
print(vis_x)
print(vis_y)

#df['tSNEx'] = vis_x
#df['tSNEy'] = vis_y
#df['variance'] = variance
#df['classif'] = classification_vector

#print(df1)
output_df = pd.concat([df,df1.transpose()])

(output_df.transpose()).to_csv(os.path.dirname(sys.argv[1])+"\\d3data_bin.tsv",sep="\t") #THIS CAN BE LOADED INTO TSNE ON D3.js