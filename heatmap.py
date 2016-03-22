__author__ = 'KATRINA'


from pandas import DataFrame
import pandas as pd
import collections
import heapq
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


input_file = "C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Output\\week2\\backSPIN\\N0_S1_R1_matrix_clustered.cef_forHeatmap"
df = DataFrame.from_csv(input_file, sep="\t")

def removeAllZeros(df):
    df_no_zeros = df.loc[(df!=0).any(axis=1)] #remove all zero-only rows
    return df_no_zeros

def removeRowsLessThanTen(df):
    df_no_under_ten = df.loc[(df.sum(axis=1) >= 10)] #remove all rows that have sum < 10
    return df_no_under_ten

def removeCols(df):
    df_no_col_under_twothousand = df.loc[:,(df.sum(axis=0) >= 1000)]
    return df_no_col_under_twothousand

df_clean = removeCols(removeRowsLessThanTen(df))
#plt.pcolor(df)
sns.heatmap(df_clean,vmin=0, vmax=600)# yticklabels=False, xticklabels=False,
#plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
#plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns)
plt.show()
