__author__ = 'KATRINA'

'''

first attempt at implementing TSNE
http://scikit-learn.org/stable/auto_examples/manifold/plot_lle_digits.html

'''

from time import time
import matplotlib.pyplot as plt
from sklearn import (manifold)
from pandas import DataFrame

df = DataFrame.from_csv("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\test_data\\CSF\\combined_matrix.txt_df_result_top20.txt", sep="\t")
assert isinstance(df, object)
X = df.transpose()
n_samples, n_features = X.shape[0],X.shape[1]
n_neighbors = 30

# t-SNE embedding of the digits dataset
print("Computing t-SNE embedding")
tsne = manifold.TSNE(n_components=2, init='pca', random_state=0)
t0 = time()
X_tsne = tsne.fit_transform(X)

# plot the result
vis_x = X_tsne[:, 0]
vis_y = X_tsne[:, 1]

plt.scatter(vis_x, vis_y)
plt.show()
