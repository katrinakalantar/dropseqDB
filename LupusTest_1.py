__author__ = 'KATRINA'

import pandas as pd
#import sklearn.feature_selection as fs
import sklearn
from sklearn import svm, cross_validation, feature_selection
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.pipeline import make_pipeline
from sklearn.pipeline import Pipeline
import matplotlib as plt
import numpy as np
import os

def normalizeRPM(df):
    total_transcripts_per_cell = []
    for i in df.columns.values:
        sum = 0
        for genes in df[i]:
            sum += genes
        total_transcripts_per_cell.append(sum)

    total_sum = np.sum(total_transcripts_per_cell)
    per_million_scaling_factor = total_sum/1000000
    new_df = df.apply(lambda x: (x/per_million_scaling_factor))
    return new_df

prefix = "C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\Lupus_data"

df1 = normalizeRPM(pd.read_csv(os.path.join(prefix,"P-19_S9_L007_R1_001_matrix.txt"),sep="\t", index_col=0))
df2 = normalizeRPM(pd.read_csv(os.path.join(prefix,"P-25_S11_L007_R1_001_matrix.txt"),sep="\t", index_col=0))
df3 = normalizeRPM(pd.read_csv(os.path.join(prefix,"P-27_S12_L007_R1_001_matrix.txt"),sep="\t", index_col=0))
df4 = normalizeRPM(pd.read_csv(os.path.join(prefix,"PF-25_S13_L007_R1_001_matrix.txt"),sep="\t", index_col=0))
df5 = normalizeRPM(pd.read_csv(os.path.join(prefix,"PW-25_S14_L007_R1_001_matrix.txt"),sep="\t", index_col=0))
df6 = normalizeRPM(pd.read_csv(os.path.join(prefix,"PW-27_S15_L007_R1_001_matrix.txt"),sep="\t", index_col=0))
df7 = normalizeRPM(pd.read_csv(os.path.join(prefix,"NID0157_dsq_1_S1_R1_matrix.txt"),sep="\t", index_col=0))
df8 = normalizeRPM(pd.read_csv(os.path.join(prefix,"NID0157_dsq_2_S4_R1_matrix.txt"),sep="\t", index_col=0))
df9 = normalizeRPM(pd.read_csv(os.path.join(prefix,"NID0194_dsq_1_S2_R1_matrix.txt"),sep="\t", index_col=0))
df10 = normalizeRPM(pd.read_csv(os.path.join(prefix,"NID0194_dsq_2_S3_R1_matrix.txt"),sep="\t", index_col=0))

result = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9,df10], axis=1)
result.to_csv(os.path.join("/data/katrina/Data/Lupus_analysis/","combined_matrix.txt"),sep="\t")


#result = pd.read_csv(os.path.join(prefix,"combined_matrix.txt"),sep="\t")
print("saved combined df file")

df = result
df_trans = result.transpose()

classification = {}
classification["SampleID"] = [0]*(len(df1.columns.values)+len(df2.columns.values)+len(df3.columns.values)) + [1]*(len(df4.columns.values)+len(df5.columns.values)+len(df6.columns.values)) +\
[2]*(len(df7.columns.values)+len(df8.columns.values)) + [3]*(len(df9.columns.values)+len(df10.columns.values))
classification["Diagnosis"] = [0]*(len(df1.columns.values)+len(df2.columns.values)+len(df3.columns.values)) + [1]*(len(df4.columns.values)+len(df5.columns.values)+len(df6.columns.values)) +\
[1]*(len(df7.columns.values)+len(df8.columns.values)) + [1]*(len(df9.columns.values)+len(df10.columns.values))

sample_classification = classification["Diagnosis"]

'''
# select k best - TODO: Look up how this algorithm works behind the scenes
X_new = fs.SelectKBest(fs.chi2, k=20).fit_transform(df_trans, sample_classification)
n = fs.SelectKBest(fs.chi2,k=20).fit(df_trans, sample_classification)

#print(X_new)
print(n)
t = n.get_support()
indices_of_interest = [i for i, x in enumerate(t) if x]
genes_of_interest = []

for i in indices_of_interest:
    genes_of_interest.append(df_trans.columns.values[i])

print(genes_of_interest)
subset_df = df[df.index.isin(genes_of_interest)]
subset_df.to_csv("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\Lupus_data\\" + "_df_result.txt",sep="\t")
'''

print("beginning feature-selection")

# Create a feature-selection transform and an instance of SVM that we
# combine together to have an full-blown estimator

transform = feature_selection.SelectPercentile(feature_selection.f_classif)

clf = Pipeline([('anova', transform), ('svc', svm.SVC(C=1.0))])

###############################################################################
# Plot the cross-validation score as a function of percentile of features
score_means = list()
score_stds = list()
percentiles = (1, 3, 6, 10, 15, 20, 30, 40, 60, 80, 100)

for percentile in percentiles:
    print("percentile: " + str(percentile))
    clf.set_params(anova__percentile=percentile)
    # Compute cross-validation score using all CPUs
    this_scores = cross_validation.cross_val_score(clf, df_trans, sample_classification, n_jobs=1)
    score_means.append(this_scores.mean())
    score_stds.append(this_scores.std())


plt.errorbar(percentiles, score_means, np.array(score_stds))

plt.title(
    'Performance of the SVM-Anova varying the percentile of features selected')
plt.xlabel('Percentile')
plt.ylabel('Prediction rate')

plt.axis('tight')
#plt.show()
plt.savefig(os.path.join("/data/katrina/Data/Lupus_analysis/","classification_plot.png"))



'''
# ANOVA SVM-C
# 1) anova filter, take 3 best ranked features
anova_filter = SelectKBest(fs.chi2, k=20)
# 2) svm
clf = svm.SVC(kernel='linear')

anova_svm = make_pipeline(anova_filter, clf)
anova_svm.fit(df_trans, sample_classification)
anova_svm.predict(df_trans)
'''