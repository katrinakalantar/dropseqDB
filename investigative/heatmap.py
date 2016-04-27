__author__ = 'KATRINA'

'''

create a heatmap from input dataframe
note: I was using this specifcially with output of BackSPIN algorithm

python heatmap.py [input_file]

'''

from pandas import DataFrame
import seaborn as sns
import matplotlib.pyplot as plt
import sys

input_file = sys.argv[1]
df = DataFrame.from_csv(input_file, sep="\t")

def removeRowsLessThanTen(df):
    df_no_under_ten = df.loc[(df.sum(axis=1) >= 10)] #remove all rows that have sum < 10
    return df_no_under_ten

def removeCols(df):
    df_no_col_under_twothousand = df.loc[:,(df.sum(axis=0) >= 1000)]
    return df_no_col_under_twothousand

df_clean = removeCols(removeRowsLessThanTen(df))
sns.heatmap(df_clean,vmin=0, vmax=600)# yticklabels=False, xticklabels=False,
#plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
#plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns)
plt.show()
