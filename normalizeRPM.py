__author__ = 'KATRINA'

'''
Normalize an input file (hard-coded)

This script contains both the original method for normalization (full file RPM), which
is actually not correct.

AND the updated normalization method "normalizeRPM_new" which normalizes per cell (column)

'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys


#normalize per cell (column) RPM
def normalizeRPM_new(df):
    total_transcripts_per_cell = 0
    for i in df.columns.values:
        sum = 0
        total_transcripts_per_cell = 0
        for genes in df[i]:
            sum += genes
        total_transcripts_per_cell=sum
        print(total_transcripts_per_cell)
        per_million_scaling_factor = float(total_transcripts_per_cell)/float(1000000)
        print(per_million_scaling_factor)
        for gene_name in df.index:
            g = float(df.ix[gene_name,i])/float(per_million_scaling_factor)
            df.set_value(gene_name, i, g)

    return df #new_df

#normalize entire file to RPM - THIS IS NOT A GOOD WAY TO NORMALIZE
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

prefix = sys.argv[1]
df1 = normalizeRPM_new((pd.read_csv(os.path.join(prefix,"NID0157_dsq_1_S1_R1_matrix.txt"),sep="\t", index_col=0)))

df1.to_csv(os.path.join(prefix,"NID0157_dsq_1_S1_R1_matrix_norm"),sep="\t")