__author__ = 'KATRINA'

'''
Generate a hierarchical clustering heatmap of the input dataframe

python clustermap_1.py [input_file]

note: in a sparse matrix this will take a long time to execute and provide
minimal useful information. Best to used this on a sub-selected matrix.

#investigative

'''

import pandas as pd
import os
import sys
import scipy
import pylab
import scipy.cluster.hierarchy as sch

input_file = sys.argv[1]
df = pd.DataFrame.from_csv((input_file),sep="\t").transpose() #because the correlation coefficient is taken for columns, we want genes on cols

# Generate features and distance matrix.
x = scipy.rand(40)
D = scipy.zeros([40,40])
for i in range(40):
    for j in range(40):
        D[i,j] = abs(x[i] - x[j])

# Compute and plot dendrogram.
fig = pylab.figure()
axdendro = fig.add_axes([0.09,0.1,0.2,0.8])
Y = sch.linkage(D, method='centroid')
Z = sch.dendrogram(Y, orientation='right')
axdendro.set_xticks([])
axdendro.set_yticks([])

# Plot distance matrix.
axmatrix = fig.add_axes([0.3,0.1,0.6,0.8])
index = Z['leaves']
D = D[index,:]
D = D[:,index]
im = axmatrix.matshow(D, aspect='auto', origin='lower')
axmatrix.set_xticks([])
axmatrix.set_yticks([])

# Plot colorbar.
axcolor = fig.add_axes([0.91,0.1,0.02,0.8])
pylab.colorbar(im, cax=axcolor)

# Display and save figure.
fig.show()
fig.savefig(os.path.join(os.path.dirname(sys.argv[1]),'dendrogram.png'))

