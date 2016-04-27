__author__ = 'KATRINA'

'''
This script was used to investigate the data associated with lupus v healthy. It was not intended to
perform a particular task or return a particular result, but rather was continually modified to
try to ways of manipulating the data

#outdated

'''

import pandas as pd
import numpy as np
import os
import json
import sys


def normalizeRPM_new(df):
    total_transcripts_per_cell = 0#[]
    for i in df.columns.values:
        sum = 0
        total_transcripts_per_cell = 0
        for genes in df[i]:
            sum += genes
        total_transcripts_per_cell=sum#.append(sum)
        #print(total_transcripts_per_cell)
        per_million_scaling_factor = float(total_transcripts_per_cell)/float(1000000)
        #print(per_million_scaling_factor)
        for gene_name in df.index:
            g = float(df.ix[gene_name,i])/float(per_million_scaling_factor)
            df.set_value(gene_name, i, g)
    return df

def normalizeRPM(df):
    total_transcripts_per_cell = []
    for i in df.columns.values:
        sum = 0
        for genes in df[i]:
            sum += genes
        total_transcripts_per_cell.append(sum)

    total_sum = np.sum(total_transcripts_per_cell)
    print("total_sum: " + str(total_sum))
    per_million_scaling_factor = float(total_sum)/float(1000000)
    print("per_million_scaling_facor: " + str(per_million_scaling_factor))
    new_df = df.apply(lambda x: (x/per_million_scaling_factor))
    return new_df


def checkNormalization(df):
    total_transcripts_per_cell = []
    for i in df.columns.values:
        sum = 0
        for genes in df[i]:
            sum += genes
        total_transcripts_per_cell.append(sum)

    total_sum = np.sum(total_transcripts_per_cell)
    print(total_sum)

def removeRowsLessThanTen(df):
    df_no_under_ten = df.loc[(df.sum(axis=1) >= 10)] #remove all rows that have sum < 10
    return df_no_under_ten

def removeCols(df):
    df_no_col_under_twothousand = df.loc[:,(df.sum(axis=0) >= 2000)]
    return df_no_col_under_twothousand

def removeZeroCols(df):
    df_no_col_under_twothousand = df.loc[:,(df.sum(axis=0) > 0)]
    return df_no_col_under_twothousand

prefix = sys.argv[1] #"C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\Lupus_data"

print("reading in files")

'''
df1 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"P-19_S9_L007_R1_001_matrix.txt"),sep="\t", index_col=0)))
df2 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"P-25_S11_L007_R1_001_matrix.txt"),sep="\t", index_col=0)))
df3 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"P-27_S12_L007_R1_001_matrix.txt"),sep="\t", index_col=0)))
df4 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"PF-25_S13_L007_R1_001_matrix.txt"),sep="\t", index_col=0)))
df5 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"PW-25_S14_L007_R1_001_matrix.txt"),sep="\t", index_col=0)))
df6 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"PW-27_S15_L007_R1_001_matrix.txt"),sep="\t", index_col=0)))
df7 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"NID0157_dsq_1_S1_R1_matrix.txt"),sep="\t", index_col=0)))
df8 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"NID0157_dsq_2_S4_R1_matrix.txt"),sep="\t", index_col=0)))
df9 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"NID0194_dsq_1_S2_R1_matrix.txt"),sep="\t", index_col=0)))
df10 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"NID0194_dsq_2_S3_R1_matrix.txt"),sep="\t", index_col=0)))
result = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9,df10], axis=1)
'''


print("df1")
df1 = removeCols(normalizeRPM_new(removeZeroCols(pd.read_csv(os.path.join(prefix,"P-19_S9_L007_R1_001_matrix.txt"),sep="\t", index_col=0))))
#checkNormalization(df1)
print("df2")
df2 = removeCols(normalizeRPM_new(removeZeroCols(pd.read_csv(os.path.join(prefix,"P-25_S11_L007_R1_001_matrix.txt"),sep="\t", index_col=0))))
#checkNormalization(df2)
print("df3")
df3 = removeCols(normalizeRPM_new(removeZeroCols(pd.read_csv(os.path.join(prefix,"P-27_S12_L007_R1_001_matrix.txt"),sep="\t", index_col=0))))
#checkNormalization(df3)




'''
print("df7")
df7 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"NID0157_dsq_1_S1_R1_matrix.txt"),sep="\t", index_col=0)))
checkNormalization(df7)
print("df8")
df8 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"NID0157_dsq_2_S4_R1_matrix.txt"),sep="\t", index_col=0)))
checkNormalization(df8)
result = pd.concat([df1,df2,df3,df7,df8], axis=1)
'''
'''
print("df9")
df9 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"NID0194_dsq_1_S2_R1_matrix.txt"),sep="\t", index_col=0)))
checkNormalization(df9)
print("df10")
df10 = normalizeRPM(removeCols(pd.read_csv(os.path.join(prefix,"NID0194_dsq_2_S3_R1_matrix.txt"),sep="\t", index_col=0)))
checkNormalization(df10)
result = pd.concat([df1,df2,df3,df9,df10], axis=1)
'''


#WAS USING THIS BLOCK
print("df4")
df4 = removeCols(normalizeRPM_new(removeZeroCols(pd.read_csv(os.path.join(prefix,"PF-25_S13_L007_R1_001_matrix.txt"),sep="\t", index_col=0))))
#checkNormalization(df4)
print("df5")
df5 = removeCols(normalizeRPM_new(removeZeroCols(pd.read_csv(os.path.join(prefix,"PW-25_S14_L007_R1_001_matrix.txt"),sep="\t", index_col=0))))
#checkNormalization(df5)
print("df6")
df6 = removeCols(normalizeRPM_new(removeZeroCols(pd.read_csv(os.path.join(prefix,"PW-27_S15_L007_R1_001_matrix.txt"),sep="\t", index_col=0))))
#checkNormalization(df6)
result = pd.concat([df1,df2,df3,df4,df5,df6], axis=1)


print("combined result dataframe has been generated")
result.to_csv(os.path.join(prefix,"combined_matrix.txt"),sep="\t")
print("saved combined df file")

df = result

#df = pd.read_csv(os.path.join(prefix,"combined_matrix.txt"),sep="\t",index_col=0)

#MODIFY DF TO REMOVE GENES WITH "OUTLIER" VARIANCE
df_trans1 = df.transpose()
print(df_trans1)
variance = df_trans1.var(axis=1) #variance in columns (genes)
sum = df_trans1.sum(axis=1) #sum of all in column - so that we can get the mean expression
print("variance array:")
print(variance.as_matrix())
print("sum array:")
print(sum.as_matrix())
var_metric = []
for i in range(len(df_trans1.columns.values)):
    #print(df1[df1.columns.values[i]])
    mean_expn = sum[i]/len(df_trans1.columns.values)
    print("mean expn:")
    print(mean_expn)
    var_metric.append((variance[i]*variance[i])/(mean_expn*mean_expn))
    print("var_matric:")
    print(var_metric)


mean = np.mean(variance)
mean_var_metric = np.mean(var_metric)
sd = np.std(variance)
sd_var_metric = np.std(var_metric)
print("var: " +str(var_metric))
print("mean: " +str(mean_var_metric))
print("std: " + str(sd_var_metric))

count_var = 0
count_sum = 0
removed_genes_var = []
col_names = df_trans1.columns.values
print("#genes before removing genes: " + str(len(df_trans1.columns.values)))
for g in col_names:
    #if (df_trans1[g].var() > 1000000): #REMOVE GENES WITH HIGH VARIANCE
    #    df_trans1.drop(g,axis=1)
    print("variance: " + str(df_trans1.as_matrix([g]).var()))
    print(df_trans1.as_matrix([g]))
    '''
    if (abs(mean - df_trans1[g].var()) > 2*sd ):
        df_trans1 = df_trans1.drop(g,axis=1)
        print("removing based on standard deviations")
        count_var +=1
        removed_genes_var.append(g)
        '''
    var = df_trans1.as_matrix([g]).var()
    sum = df_trans1.as_matrix([g]).sum()
    mean_expn = sum/len(df_trans1.as_matrix([g]))
    var_metric = (variance[i]*variance[i])/(mean_expn*mean_expn)

    #comput the most common variance metric outlier
    if (abs(mean_var_metric - var_metric) > 2*sd_var_metric):
        df_trans1 = df_trans1.drop(g,axis=1)
        print("removing based on standard deviations")
        count_var +=1
        removed_genes_var.append(g)
    elif (df_trans1[g].sum() < 2000): #REMOVE GENES WITH TOTAL COUNTS < 2000
        print(df_trans1[g].sum())
        df_trans1 = df_trans1.drop(g,axis=1)
        count_sum +=1

print("#genes after removing genes: " + str(len(df_trans1.columns.values)))
print("total genes removed due to variance: " + str(count_var))
print(removed_genes_var)
print("total genes removed due to sum: " + str(count_sum))
df = df_trans1.transpose()

print("cleaned combined result dataframe has been generated")
df.to_csv(os.path.join(prefix,"combined_matrix_limitedVar.txt"),sep="\t")
print("saved cleaned combined df file")




df_trans = df.transpose()

classification = {}
'''
classification["SampleID"] = [0]*(len(df1.columns.values)+len(df2.columns.values)+len(df3.columns.values)) + [1]*(len(df4.columns.values)+len(df5.columns.values)+len(df6.columns.values)) +\
[2]*(len(df7.columns.values)+len(df8.columns.values)) + [3]*(len(df9.columns.values)+len(df10.columns.values))
classification["Diagnosis"] = [0]*(len(df1.columns.values)+len(df2.columns.values)+len(df3.columns.values)) + [1]*(len(df4.columns.values)+len(df5.columns.values)+len(df6.columns.values)) +\
[1]*(len(df7.columns.values)+len(df8.columns.values)) + [1]*(len(df9.columns.values)+len(df10.columns.values))
'''
#classification["Diagnosis"] = [0]*(len(df1.columns.values)+len(df2.columns.values)+len(df3.columns.values)) + [1]*(len(df7.columns.values)+len(df8.columns.values))
#classification["Diagnosis"] = [0]*(len(df1.columns.values)+len(df2.columns.values)+len(df3.columns.values)) + [1]*(len(df9.columns.values)+len(df10.columns.values))
classification["Diagnosis"] = [0]*(len(df1.columns.values)+len(df2.columns.values)+len(df3.columns.values)) + [1]*(len(df4.columns.values)+len(df5.columns.values)+len(df6.columns.values))


sample_classification = classification["Diagnosis"]

json.dump(classification, open(os.path.join(prefix, "classification.json"),'w') )
