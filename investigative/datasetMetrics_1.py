__author__ = 'KATRINA'

'''

generate a data frame with statistics on expression for each gene
options:
1. sort on total expression summed over the whole sample (all cells)
2. sort on mean expression for only cells expressing that gene

> datasetMetrics_1.py [input file]

#investigative

'''

from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import sys

#import data data
input_file = sys.argv[1] #"C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\hannas_dropseq\\N5_S7_R1_matrix.txt"
df = DataFrame.from_csv(input_file, sep="\t",index_col=0)

sampleName = input_file.split("\\")[-1].split('.')[0]

setA = []
setB = []

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

plt.hist(genes_represented_by_each_cell, color="#ffcc00",log=True)
plt.title("Number of Genes Represented per Cell")
plt.xlabel("Number of Genes")
plt.ylabel("Number of Cells")
plt.savefig(input_file.split('.')[0]+"_img1.png")

plt.hist(total_transcripts_per_cell, color="#19194d", bins=100, log=True)
plt.title("Total Transcript Count per Cell")
plt.xlabel("Total Transcript Count")
plt.ylabel("Number of Cells")
plt.savefig(input_file.split('.')[0]+"_img2.png")

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

#TODO: CREATE PLOTS FOR EACH OF THE ARRAYS ABOVE ^^^



#Calculate the variance-squared over the mean for all genes and plot
array_mean = []
array_cv = []
array_div = []
labels = []
for i in df.index: #for every row (gene)
    try:
        gene_row = df.loc[i]
        labels.append(i)
        expressed_gene_row = [x for x in gene_row if x > 0]
        mean_expr_in_expr_cells = np.mean(expressed_gene_row)
        array_mean.append(mean_expr_in_expr_cells)
        cv = scipy.stats.variation(expressed_gene_row)
        array_cv.append(cv)
        array_div.append(float(cv/mean_expr_in_expr_cells))
    except:
        "blah blah"
fig, ax = plt.subplots()
ax.scatter(array_mean, array_cv)
for i, txt in enumerate(labels):
    ax.annotate(txt, (array_mean[i],array_cv[i]))
plt.title("CV v Mean distribution")
plt.savefig(input_file+"img3.png")

'''



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



