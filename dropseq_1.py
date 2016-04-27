__author__ = 'KATRINA'

'''
Pipeline-style script containing many functions to remove low expression rows/columns (genes/cells)
and remove high variance data (genes/cells).

Saves the final dataframe as [input_file]_filtere.txt

python dropseq_1.py [input_file]

note: this script has hard-coded threshold for removal of data. It is advised that you ovserve the distributions
within a particular sample using investigateGeneVariance.py, plotGeneExprDistribution.py etc prior to running
this script and modify the parameters accordingly.

'''


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from datetime import datetime

startTime = datetime.now()
input_file = sys.argv[1]

def normalizePerCell(df):
    print("normalizing input data per cell")
    total_transcripts_per_cell = 0
    for i in df.columns.values:
        sum = 0
        total_transcripts_per_cell = 0
        for genes in df[i]:
            sum += genes
        total_transcripts_per_cell=sum
        per_million_scaling_factor = float(total_transcripts_per_cell)/float(1000000)
        for gene_name in df.index:
            g = float(df.ix[gene_name,i])/float(per_million_scaling_factor)
            df.set_value(gene_name, i, g)
    print("completed normalization  -  " + str(datetime.now() - startTime))
    print(df.shape)
    return df

def removeRows(df, threshold):
    df_no_under_ten = df.loc[(df.sum(axis=1) >= threshold)] #remove all rows that have sum < 10
    print("completed removal of genes(rows): (threshold = "+ str(threshold) +")  -  " + str(datetime.now() - startTime))
    print(df_no_under_ten.shape)
    return df_no_under_ten

def removeCols(df, threshold):
    df_no_col_under_twothousand = df.loc[:,(df.sum(axis=0) >= threshold)]
    print("completed removal of cells(columns): (threshold = "+ str(threshold) +")  -  " + str(datetime.now() - startTime))
    print(df_no_col_under_twothousand.shape)
    return df_no_col_under_twothousand


def removeZeroRows(df):
    df_no_row_under_zero = df.loc[(df.sum(axis=1) > 0)]
    print("completed removal of genes(rows): (threshold = 0)  -  " + str(datetime.now() - startTime))
    print(df_no_row_under_zero.shape)
    return df_no_row_under_zero

def removeZeroCols(df):
    df_no_col_under_zero = df.loc[:,(df.sum(axis=0) > 0)]
    print("completed removal of cells(columns): (threshold = 0)  -  " + str(datetime.now() - startTime))
    print(df_no_col_under_zero.shape)
    return df_no_col_under_zero

def removeLowExpn(df, value, threshold):
    if value == 'Gene':
        return (removeRows(df,threshold))
    elif value == 'Cell':
        return (removeCols(df,threshold))
    else:
        print("ERROR: must set value to Gene or Cell")


def viewCellExpnDistribution(df):
    genes_represented_by_each_cell = []
    total_transcripts_per_cell = []
    for i in df.columns.values:
        count = 0
        sum = 0
        for genes in df[i]:
            if genes > 0:
                count +=1
                sum +=genes
        genes_represented_by_each_cell.append(count)
        total_transcripts_per_cell.append(sum)

    plt.hist(genes_represented_by_each_cell, color="#ffcc00",alpha = 0.75)
    plt.hist(total_transcripts_per_cell, color="#19194d", bins=100,alpha = 0.75)
    plt.title("Distribution of Transcript Count per Cell")
    plt.savefig(input_file.split('.')[0]+"_countsPerCellDist.png")
    print("saved figure for Distribution of Transcript Count per Cell  -  " + str(datetime.now() - startTime))
    plt.close()

def viewGeneExpnDistribution(df):
    #distribution of reads/gene
    cells_containing_each_gene = []
    total_transcripts_per_gene = []
    variance_trancripts_per_gene = []
    for i in df.index:
        count = 0
        sum = 0
        variance = []
        for cell in df.loc[[i]]:
            x = df.loc[i,cell]
            if x > 0:
                count +=1
                sum += x
                variance.append(x)
        cells_containing_each_gene.append(count)
        total_transcripts_per_gene.append(sum)
        variance_trancripts_per_gene.append(np.var(variance))

    plt.hist(cells_containing_each_gene, color="#2bf60b")
    plt.hist(total_transcripts_per_gene, color="#19194d", bins=100)
    plt.title("Distribution of Transcript Count per Gene")
    plt.savefig(input_file.split('.')[0]+"_countsPerGeneDist.png",alpha = 0.75)
    print("saved figure for Distribution of Transcript Count per Gene  -  " + str(datetime.now() - startTime))
    plt.close()

def removeHighVariance(df): ##THIS IS SUUUUPER SLOW
    df_trans1 = df.transpose()
    print(df_trans1.shape)
    #NOTE: THESE SHOULD HAVE BEEN 1 I THOUGHT TO GET COLUMN SUMS...BUT THEY ARE CAUSING AN ERROR THAT WAY
    variance = df_trans1.var(axis=0) #variance in columns (genes)
    sum = df_trans1.sum(axis=0) #sum of all in column - so that we can get the mean expression
    var_metric = []
    print(len(df_trans1.columns.values))
    print(sum.shape)
    for i in range(len(df_trans1.columns.values)):
        print(i)
        mean_expn = sum[i]/len(df_trans1.columns.values)
        var_metric.append((variance[i])/(mean_expn*mean_expn))

    mean_var_metric = np.mean(var_metric)
    sd_var_metric = np.std(var_metric)
    print("variance metric: " +str(var_metric))
    print("mean of variance metric: " +str(mean_var_metric))
    print("std of variance metric: " + str(sd_var_metric))

    count_var = 0
    count_sum = 0
    removed_genes_var = []
    col_names = df_trans1.columns.values
    print("#genes before removing genes: " + str(len(df_trans1.columns.values)))
    for g in col_names:
        var = df_trans1.as_matrix([g]).var()
        sum = df_trans1.as_matrix([g]).sum()
        mean_expn = sum/len(df_trans1.as_matrix([g]))
        var_metric = (variance[i]/(mean_expn*mean_expn))
        #compute the most common variance metric outlier
        if (abs(mean_var_metric - var_metric) > 2*sd_var_metric):
            df_trans1 = df_trans1.drop(g,axis=1)
            print("removing based on standard deviations")
            count_var +=1
            removed_genes_var.append(g)

    print("#genes after removing genes: " + str(len(df_trans1.columns.values)))
    print("total genes removed due to variance: " + str(count_var))
    print(removed_genes_var)
    print("total genes removed due to sum: " + str(count_sum))
    df = df_trans1.transpose()
    print("completed removal of genes with high variance  -  " + str(datetime.now() - startTime))
    print(df.shape)
    return df

df_original = pd.read_csv(input_file,sep="\t", index_col=0)
print("read in original data to data frame  -  " + str(datetime.now() - startTime))
print(df_original.shape)
#df = removeHighVariance(removeRows(removeCols(normalizePerCell(removeZeroRows(removeZeroCols(df_original))),10000),5000))

df_rc = removeCols(df_original, 2000)
df = normalizePerCell(df_rc)
df.to_csv(input_file.split('.')[0]+"_normalized.txt",sep="\t")

df_final = removeRows(removeCols(df_original,10000),5000)
#viewCellExpnDistribution(df) #this is dumb because they are all around 0 or 1 million after normalization
viewGeneExpnDistribution(df)

df_final.to_csv(input_file.split('.')[0]+"_filtered.txt",sep="\t")
