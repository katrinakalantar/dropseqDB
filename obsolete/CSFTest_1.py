__author__ = 'KATRINA'

'''
Playing around with CSF data -

loads all data into a dataframe and creates the classification vector
the combined data is saved in a CSF_combined_matrix.txt file and the classification vector
is saved in classification.json

'''

import pandas as pd
import numpy as np
import os
import json
import sys

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

def removeRowsLessThanTen(df):
    df_no_under_ten = df.loc[(df.sum(axis=1) >= 10)] #remove all rows that have sum < 10
    return df_no_under_ten

def removeCols(df):
    df_no_col_under_twothousand = df.loc[:,(df.sum(axis=0) >= 1000)]
    return df_no_col_under_twothousand


prefix = sys.argv[1] #"C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\Lupus_data"

print("reading in files")
df1 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"Cyro_UnS_0_1_S20_R1_matrix.txt"),sep="\t", index_col=0)))
df2 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"Cyro_UnS_0_2_S21_R1_matrix.txt"),sep="\t", index_col=0)))
df3 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"N0_S1_R1_matrix.txt"),sep="\t", index_col=0)))
df4 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"N1_S2_R1_matrix.txt"),sep="\t", index_col=0)))
df5 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"N2_S3_R1_matrix.txt"),sep="\t", index_col=0)))
df6 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"N3_S5_R1_matrix.txt"),sep="\t", index_col=0)))
df7 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"N5_S7_R1_matrix.txt"),sep="\t", index_col=0)))

result = pd.concat([df1,df2,df3,df4,df5,df6,df7], axis=1)
print("combined result dataframe has been generated")
print(prefix)
result.to_csv(os.path.join(prefix,"CSF_combined_matrix.txt"),sep="\t")

#result = pd.read_csv(os.path.join(prefix,"combined_matrix.txt"),sep="\t")
print("saved combined df file")

df = result
df_trans = result.transpose()

classification = {}
classification["SampleID"] = [0]*(len(df1.columns.values)+len(df2.columns.values)) + [1]*(len(df3.columns.values))+ \
                             [2]*(len(df4.columns.values))+ [3]*(len(df5.columns.values))+ [4]*(len(df6.columns.values))+ \
                             [5]*(len(df7.columns.values))
#sample_classification = classification["SampleID"]

json.dump(classification, open(os.path.join(prefix, "classification.json"),'w') )
