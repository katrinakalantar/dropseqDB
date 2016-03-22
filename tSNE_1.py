'''
http://scikit-learn.org/stable/auto_examples/manifold/plot_lle_digits.html

'''

__author__ = 'KATRINA'


from time import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import offsetbox
from sklearn import (manifold, datasets, decomposition, ensemble, random_projection)
from pandas import DataFrame

#digits = datasets.load_digits(n_class=6)
#X = digits.data
#y = digits.target
#n_samples, n_features = X.shape
#n_neighbors = 30

df = DataFrame.from_csv("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\test_data\\CSF\\combined_matrix.txt_df_result_top20.txt", sep="\t")
assert isinstance(df, object)
X = df.transpose()
n_samples, n_features = X.shape[0],X.shape[1]
n_neighbors = 30

'''
#----------------------------------------------------------------------
# Scale and visualize the embedding vectors
def plot_embedding(X, title=None):
    x_min, x_max = np.min(X, 0), np.max(X, 0)
    X = (X - x_min) / (x_max - x_min)

    plt.figure()
    ax = plt.subplot(111)

    plt.xticks([]), plt.yticks([])
    if title is not None:
        plt.title(title)
'''


#----------------------------------------------------------------------
# t-SNE embedding of the digits dataset
print("Computing t-SNE embedding")
tsne = manifold.TSNE(n_components=2, init='pca', random_state=0)
t0 = time()
X_tsne = tsne.fit_transform(X)

'''
plot_embedding(X_tsne,
               "t-SNE embedding of the digits (time %.2fs)" %
               (time() - t0))
'''

# plot the result
vis_x = X_tsne[:, 0]
vis_y = X_tsne[:, 1]

plt.scatter(vis_x, vis_y)#, cmap=plt.cm.get_cmap("jet", 10))
#plt.colorbar(ticks=range(10))
#plt.clim(-0.5, 9.5)
plt.show()
