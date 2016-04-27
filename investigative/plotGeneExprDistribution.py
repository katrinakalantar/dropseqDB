__author__ = 'KATRINA'

'''

plot the distributions of transcript counts across
1. genes
2. cells
3. specific gene across all cells

plotGeneExprDistribution.py [input_file]

#investigative

'''


from pandas import DataFrame
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import sys

input_file = sys.argv[1]
df = DataFrame.from_csv(input_file, sep="\t")
#df = DataFrame.from_csv("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\Lupus_data\\P-19_S9_L007_R1_001_matrix.txt", sep="\t")


def removeAllZeros(df):
    df_no_zeros = df.loc[(df!=0).any(axis=1)] #remove all zero-only rows
    return df_no_zeros

def removeRowsLessThanTen(df):
    df_no_under_ten = df.loc[(df.sum(axis=1) >= 10)] #remove all rows that have sum < 10
    return df_no_under_ten

def removeCols(df):
    df_no_col_under_twothousand = df.loc[:,(df.sum(axis=0) >= 2000)]
    return df_no_col_under_twothousand

#plot histograms...
def plotDistribution(df):

    #for transcript count across genes
    sum_of_rows = df.sum(axis=1) #sum over all rows in data frame
    plt.hist(sum_of_rows,bins=400,color="#00cc00", alpha=.5, log=True)
    plt.xlim((0,20000))
    plt.title('Distribution of Transcript Count Across Genes')
    plt.xlabel('total raw transcript count per gene')
    plt.ylabel('number of genes')
    plt.show()

    #for transcript count across cells,
    sum_of_cols = df.sum(axis=0) #sum over all columns (cells) in data frame
    plt.hist(sum_of_cols,bins=100,color="#009999", alpha=.5, log=True)
    plt.title('Distribution of Transcript Count per Cell')
    plt.xlabel('total raw transcript count within cell')
    plt.ylabel('number of cells')
    plt.show()

    #for counts across cells for an individual (randomly selected gene)
    geneA = df.loc['BLNK']
    geneB = df.loc['NCF2']
    geneC = df.loc['BCL11A']
    geneD = df.loc['VCAN']
    geneE = df.loc['NBPF11']
    geneF = df.loc['SDF4']

    #individual gene plots
    '''
    plt.hist(geneA, bins=25, alpha=.5)
    plt.title('BLNK')
    plt.xlabel('normalized transcript count per cell')
    plt.ylabel('number of cells')
    plt.show()
    plt.hist(geneB, bins=25, alpha=.5)
    plt.title('NCF2')
    plt.xlabel('normalized transcript count per cell')
    plt.ylabel('number of cells')
    plt.show()
    plt.hist(geneC, bins=25, alpha=.5)
    plt.title('BCL11A')
    plt.xlabel('normalized transcript count per cell')
    plt.ylabel('number of cells')
    plt.show()
    plt.hist(geneD, bins=25, alpha=.5)
    plt.title('VCAN')
    plt.xlabel('normalized transcript count per cell')
    plt.ylabel('number of cells')
    plt.show()
    '''

    #combined gene histogram for 3 genes
    plt.hist(geneE, bins=100, alpha=.3, color="#ff0000",log=True)
    plt.hist(geneF, bins=8, alpha=.5, color="#0066ff",log=True)
    plt.hist(geneD, bins=100, alpha=.6, color="#ffff00",log=True)
    plt.xlim((0,2400))
    plt.xlabel('normalized transcript count per cell', fontsize=16)
    plt.ylabel('number of cells',fontsize=16)

    yellow_patch = mpatches.Patch(color='#ffff00', label='VCAN')
    red_patch = mpatches.Patch(color='#ff0000', label='NBPF11')
    blue_patch = mpatches.Patch(color='#0066ff',label='SDF4')
    plt.legend(handles=[yellow_patch,red_patch,blue_patch])
    plt.show()


#df_clean = removeCols(removeRowsLessThanTen(df))
plotDistribution(df)
