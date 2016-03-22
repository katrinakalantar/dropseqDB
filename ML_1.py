__author__ = 'KATRINA'

'''

investigating Machine Learning (ML) methods for determining the most informative genes for separating two datasets
python ML_1.py inputDF classificationFile [-norm]

'''

from time import time
from pandas import DataFrame
import pandas as pd
import collections
import heapq
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import sklearn.feature_selection as fs
import json
import sys
import mygene


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

input_file = sys.argv[1] #"C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\test_data\\ML_test\\test_out_RPM.tsv"
df = pd.read_csv(input_file,sep="\t", index_col=0)
if '-norm' in sys.argv:
    df_norm = normalizeRPM(df)
    df = df_norm


print(len(df.columns.values))
print(len(df.index))
df_trans = df.transpose()

sample_classification = json.load(open(sys.argv[2]))["SampleID"]#"C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\test_data\\ML_test\\test_outclassification.txt"))["SampleID"] #
print(sample_classification)
print(len(sample_classification))

# select k best - TODO: Look up how this algorithm works behind the scenes
X_new = fs.SelectKBest(fs.chi2, k=1000).fit_transform(df_trans, sample_classification)
n = fs.SelectKBest(fs.chi2,k=1000).fit(df_trans, sample_classification)

print(n)
t = n.get_support()
indices_of_interest = [i for i, x in enumerate(t) if x]
genes_of_interest = []

for i in indices_of_interest:
    genes_of_interest.append(df_trans.columns.values[i])

print(genes_of_interest)
subset_df = df[df.index.isin(genes_of_interest)]
subset_df.to_csv(sys.argv[1]+"_df_result.txt",sep="\t")




##TESTING - GET PROTEIN INFORMATION FROM DATABASE, USE THIS WITH CYTOSCAPE
mg = mygene.MyGeneInfo()
a = mg.querymany(genes_of_interest, scopes="symbol", fields = ["uniprot", "ensembl.gene", "reporter"], species="human", as_dataframe=True)
#print(a)
#print(a['ensembl.gene'])
a['ensembl.gene'].to_csv(sys.argv[1]+"_DFgeneSymbols")

'''
for i in genes_of_interest:
    print('symbol:'+i)
    mg.query('symbol:'+i, species='human')
'''