__author__ = 'KATRINA'


'''

generate a data frame with statistics on expression for each gene
options:
1. sort on total expression summed over the whole sample (all cells)
2. sort on mean expression for only cells expressing that gene

#investigative

'''

__author__ = 'KATRINA'

from time import time
from pandas import DataFrame
import pandas as pd
import collections
import heapq
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

#real data
input_file = "C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\hannas_dropseq\\N0_S1_R1_matrix.txt"
df = DataFrame.from_csv(input_file, sep="\t")

sampleName = input_file.split("\\")[-1].split('.')[0]

#col_names = list(df.columns.values)
#row_names = df.index

def removeAllZeros(df):
    df_no_zeros = df.loc[(df!=0).any(axis=1)] #remove all zero-only rows
    return df_no_zeros

def removeRowsLessThanTen(df):
    df_no_under_ten = df.loc[(df.sum(axis=1) >= 10)] #remove all rows that have sum < 10
    return df_no_under_ten

def removeCols(df):
    df_no_col_under_twothousand = df.loc[:,(df.sum(axis=0) >= 500)]
    return df_no_col_under_twothousand


setA = []
setB = []


#take only subset of data that meets initial filtering conditions
df_clean = removeCols(removeRowsLessThanTen(df))
print("original file #cells: " + str(len(df.columns.values)))   #print statistics on what was filtered out
print("original file #genes: " + str(len(df.index)))
print("filtered file #cells: " + str(len(df_clean.columns.values)))
print("filtered file #genes: " + str(len(df_clean.index)))

gene_sums_across_cells = df_clean.sum(axis=1)
sns.distplot(gene_sums_across_cells)   #plot distribution of counts
sns.plt.show()

plt.hist(gene_sums_across_cells,bins=250)
plt.axis([0, 2000, 0, 6000])
plt.show()

row_names = df_clean.index
totalCounts = np.sum(gene_sums_across_cells)

gene_percents = [x/totalCounts for x in gene_sums_across_cells]
print(gene_percents)
sns.distplot(gene_percents)    #plot distribution of percents (should be the same as the counts)
sns.plt.show()

plt.hist(gene_percents,bins=50)
plt.axis([0, .1, 0, 100])
plt.show()

gene_counts = collections.defaultdict(list) #create dictionary of gene counts : gene names
for sum,name in zip(gene_sums_across_cells, row_names):
    gene_counts[sum].append(name)
print(gene_counts)


largest_gene_counts = heapq.nlargest(30,gene_counts)
totalNumCells = len(df_clean.columns.values)

output_DF1 = DataFrame(columns=['Gene Name', 'Count', 'Fraction_of_Total', 'Fraction_of_cells_expressing_gene', 'Mean_expn_in_expressing_cells', 'Std_expn_in_expression_cells'])
output_DF2 = DataFrame(columns=['Gene Name', 'Count', 'Fraction_of_Total', 'Fraction_of_cells_expressing_gene', 'Mean_expn_in_expressing_cells', 'Std_expn_in_expression_cells'])

#SORTING ON TOTAL EXPRESSION FOR THE ENTIRE SAMPLE
for i in largest_gene_counts: #loop through integer values of largest counts
    PofT = i / totalCounts
    for g in gene_counts[i]:
        setA.append(g)
        gene_row = df_clean.loc[g]
        expressed_gene_row = [x for x in gene_row if x > 0]
        fraction_cells_expr_gene = len(expressed_gene_row)/totalNumCells
        mean_expr_in_expr_cells = np.mean(expressed_gene_row)
        std_expr_in_expr_cells = np.std(expressed_gene_row)
        #print the final outcome
        output_DF1.loc[len(output_DF1)] = [g,str(i),str(np.round(PofT,decimals=3)),str(np.round(fraction_cells_expr_gene,decimals=3)),str(np.round(mean_expr_in_expr_cells,decimals=3)),str(np.round(std_expr_in_expr_cells,decimals=3))]


#SORTING ON MEAN EXPRESSION WITHIN EXPRESSING CELLS
for i in df_clean.index: #for every row (gene)
    gene_row = df_clean.loc[i]
    expressed_gene_row = [x for x in gene_row if x > 0]
    fraction_cells_expr_gene = len(expressed_gene_row)/totalNumCells
    mean_expr_in_expr_cells = np.mean(expressed_gene_row)
    std_expr_in_expr_cells = np.std(expressed_gene_row)
    gene_count = np.sum(gene_row)

    output_DF2.loc[len(output_DF2)]=[i,gene_count,np.round(gene_count/totalCounts,decimals=3),np.round(fraction_cells_expr_gene,decimals=3),np.round(mean_expr_in_expr_cells,decimals=3),np.round(std_expr_in_expr_cells,decimals=3)]

new_DF1 = output_DF1.sort_index(by='Count',ascending=0)
new_DF1.to_csv("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Output\\"+sampleName+"topGenes_1.output.tsv",sep='\t')

new_DF2 = output_DF2.sort_index(by='Mean_expn_in_expressing_cells',ascending=0)[0:30]
new_DF2.to_csv("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Output\\"+sampleName + "topGenes_2.output.tsv",sep='\t')


#get overlapping set from top 50 of each
a = new_DF1['Gene Name'][0:50]
for s in a:
    setA.append(s)
b = new_DF2['Gene Name'][0:50]
for s in b:
    setB.append(s)

set_A = set(setA)
set_B = set(setB)

final_gene_set = set_A.intersection(set_B)
print(final_gene_set)


'''
TODO: This part doesn't work - need a way to access the data frame values for correct genes
output_final = DataFrame(columns=['Gene Name', 'Count', 'Fraction_of_Total', 'Fraction_of_cells_expressing_gene', 'Mean_expn_in_expressing_cells', 'Std_expn_in_expression_cells'])
for i in final_gene_set:
    output_final.loc[len(output_final)]=[i,gene_count,gene_count/totalCounts,fraction_cells_expr_gene,mean_expr_in_expr_cells,std_expr_in_expr_cells]
'''



