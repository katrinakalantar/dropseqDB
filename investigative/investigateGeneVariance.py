__author__ = 'KATRINA'

'''
Plot the distribution of variance across genes and cells.

investigateGeneVariance.py [file of interest]

#investigative

'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import math

#remove outliers from the data frame
def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]

#prefix = sys.argv[1]
#df1 = pd.read_csv(os.path.join(prefix,"NID0157_dsq_1_S1_R1_matrix.txt"),sep="\t", index_col=0)

filename = sys.argv[1]
df1 = pd.read_csv(filename, sep="\t", index_col=0)

#get the variance in rows (gene expression), plot it
variance = df1.var(axis=1) #variance in rows
sum = df1.sum(axis=1)
a = []
var_metric = []
for i in range(len(variance)):
    a.append(variance[i])
    mean_expn = sum[i]/len(df1.columns.values)
    x = (variance[i]*variance[i])/(mean_expn*mean_expn)
    if not math.isnan(x):
        var_metric.append((variance[i]*variance[i])/(mean_expn*mean_expn))

#plot the variance in gene expression prior to removal of outliers
#plt.hist(var_metric, color="#ff0066",bins=np.logspace(0.1, 8, 50),alpha = 0.75)
#plt.gca().set_xscale("log")
#plt.show()
#plt.close()

x = reject_outliers(np.asarray(a))
plt.hist(x,color="#ff3300", bins=np.logspace(0.1, 8, 50),alpha=0.5)#bins=50, alpha=0.5)
plt.gca().set_xscale("log")
plt.title("distribution of variance in genes expression (rows) - outliers removed")
plt.show()
plt.hist(variance[variance!=0])
plt.close()


#get the variance in cols (cells), plot it
a = []
cell_to_cell_var = df1.var(axis=0) #variance in cols (cells)
for i in cell_to_cell_var:
    a.append(i)
def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]
x = reject_outliers(np.asarray(a))
plt.hist(x,color="#0099ff", bins=50, alpha=0.5)
plt.title("distribution of variance in cell expression (cols) - outliers removed")
plt.show()
