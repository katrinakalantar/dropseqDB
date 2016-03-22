'''

create a histogram of total counts per gene
#investigative

'''

__author__ = 'KATRINA'

from time import time
from pandas import DataFrame
import seaborn as sns
import numpy as np


import matplotlib.pyplot as plt
from matplotlib import offsetbox
from sklearn import (manifold, datasets, decomposition, ensemble, random_projection)




#testing
#df = DataFrame.from_csv("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\test_data\\N0_S1_R1_matrixTi.txt", sep="\t")
df = DataFrame.from_csv("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\test_data\\N1_S2_R1_matrixSm.txt", sep="\t")

#real data
#df = DataFrame.from_csv("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\hannas_dropseq\\N0_S1_R1_matrix.txt", sep="\t")

#col_names = list(df.columns.values)
#row_names = df.index

def removeAllZeros(df):
    df_no_zeros = df.loc[(df!=0).any(axis=1)] #remove all zero-only rows
    return df_no_zeros

def removeRowsLessThanTen(df):
    df_no_under_ten = df.loc[(df.sum(axis=1) >= 10)] #remove all rows that have sum < 10
    return df_no_under_ten

def removeCols(df):
    df_no_col_under_twothousand = df.loc[:,(df.sum(axis=0) >= 2000)]
    return df_no_col_under_twothousand

#plot histogram and kde plot
def plotDistribution(df):
    sum_of_rows = df.sum(axis=1) #sum over all rows in data frame
    sns.distplot(sum_of_rows, bins=100)
    sns.plt.show()
    sns.kdeplot(sum_of_rows)
    sns.plt.show()



def scatter(x, colors):
    # We choose a color palette with seaborn.
    palette = np.array(sns.color_palette("hls", 10))

    # We create a scatter plot.
    f = plt.figure(figsize=(8, 8))
    ax = plt.subplot(aspect='equal')
    sc = ax.scatter(x[:,0], x[:,1], lw=0, s=40)#, c=palette[colors.astype(np.int)])
    plt.xlim(-25, 25)
    plt.ylim(-25, 25)
    ax.axis('off')
    ax.axis('tight')


    return f, ax, sc


def plotTSNE(df):
    X = df
    n_samples, n_features = X.shape[0],X.shape[1]
    n_neighbors = 30
    print("Computing t-SNE embedding")
    tsne = manifold.TSNE(n_components=2, init='pca', random_state=0)
    t0 = time()
    X_tsne = tsne.fit_transform(X)
    print("finished tsne.fit_transform()")


    # plot the result
    #vis_x = X_tsne[:, 0]
    #vis_y = X_tsne[:, 1]
    print('blah')

    '''
    plt.scatter(vis_x, vis_y, cmap=plt.cm.get_cmap("jet", 10))
    plt.colorbar(ticks=range(10))
    plt.clim(-0.5, 9.5)
    plt.show()
    '''
    #y = np.hstack([df.target[df.target==i] for i in range(10)])
    y=0
    scatter(X_tsne,y)
    plt.show()


df_clean = removeCols(removeRowsLessThanTen(df))
#plotDistribution(df_clean)
plotTSNE(df_clean)