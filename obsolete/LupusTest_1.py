__author__ = 'KATRINA'

'''
This script was used to investigate the data associated with lupus v healthy. It was not intended to
perform a particular task or return a particular result, but rather was continually modified to
try to ways of manipulating the data

#outdated

'''

import pandas as pd
from sklearn import svm, cross_validation, feature_selection
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
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
df1 = normalizeRPM_new(removeCols(pd.read_csv(os.path.join(prefix,"P-19_S9_L007_R1_001_matrix.txt"),sep="\t", index_col=0)))
#checkNormalization(df1)
print("df2")
df2 = normalizeRPM_new(removeCols(pd.read_csv(os.path.join(prefix,"P-25_S11_L007_R1_001_matrix.txt"),sep="\t", index_col=0)))
#checkNormalization(df2)
print("df3")
df3 = normalizeRPM_new(removeCols(pd.read_csv(os.path.join(prefix,"P-27_S12_L007_R1_001_matrix.txt"),sep="\t", index_col=0)))
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

print("df4")
df4 = normalizeRPM_new(removeCols(pd.read_csv(os.path.join(prefix,"PF-25_S13_L007_R1_001_matrix.txt"),sep="\t", index_col=0)))
#checkNormalization(df4)
print("df5")
df5 = normalizeRPM_new(removeCols(pd.read_csv(os.path.join(prefix,"PW-25_S14_L007_R1_001_matrix.txt"),sep="\t", index_col=0)))
#checkNormalization(df5)
print("df6")
df6 = normalizeRPM_new(removeCols(pd.read_csv(os.path.join(prefix,"PW-27_S15_L007_R1_001_matrix.txt"),sep="\t", index_col=0)))
#checkNormalization(df6)
result = pd.concat([df1,df2,df3,df4,df5,df6], axis=1)


print("combined result dataframe has been generated")
result.to_csv(os.path.join(prefix,"combined_matrix.txt"),sep="\t")
print("saved combined df file")

df = result

#MODIFY DF TO REMOVE GENES WITH "OUTLIER" VARIANCE #TODO: this is hard-coded to 1million, change that.
df_trans1 = df.transpose()

variance = df_trans1.var(axis=0) #variance in columns (genes)
mean = np.mean(variance)
sd = np.std(variance)

print("#genes before removing genes: " + str(len(df_trans1.columns.values)))
for g in df_trans1.columns.values:
   #if (df_trans1[g].var() > 1000000): #REMOVE GENES WITH HIGH VARIANCE
    #    df_trans1.drop(g,axis=1)
    if (abs(mean - df_trans1[g].var()) > 2*sd ):
        df_trans1.drop(g,axis=1)
        print("removing based on standard deviations")
    elif (df_trans1[g].sum() < 10): #REMOVE GENES WITH TOTAL COUNTS < 10
        df_trans1.drop(g,axis=1)
        print("removing based on sum of rows")
print("#genes after removing genes: " + str(len(df_trans1.columns.values)))

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
'''
print("selecting k best")

# select k best - TODO: Look up how this algorithm works behind the scenes
X_new = feature_selection.SelectKBest(feature_selection.chi2, k=100).fit_transform(df_trans, sample_classification)
n = feature_selection.SelectKBest(feature_selection.chi2,k=100).fit(df_trans, sample_classification)

#print(X_new)
print(n)
t = n.get_support()
indices_of_interest = [i for i, x in enumerate(t) if x]
print(indices_of_interest)
genes_of_interest = []

print("obtaining genes of interest")
for i in indices_of_interest:
    genes_of_interest.append(df_trans.columns.values[i])

print(genes_of_interest)
goi_file = open(os.path.join(prefix,"kbest_genes_of_interest.txt"),'w')
goi_file.write(str(genes_of_interest))
subset_df = df[df.index.isin(genes_of_interest)]
print(len(subset_df.columns.values))
print(len(subset_df.index))
print("writing subset dataframe to file")
subset_df.to_csv("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\Lupus_data\\" + "_df_result.tsv",sep="\t")
'''

print("beginning feature-selection")

# Create a feature-selection transform and an instance of SVM that we
# combine together to have an full-blown estimator

transform = feature_selection.SelectPercentile(feature_selection.f_classif)
anova_filter = feature_selection.SelectKBest(feature_selection.chi2)

clf = Pipeline([('anova', transform), ('svc', svm.SVC(C=1.0))])

clf_kbest = svm.SVC(kernel='linear')
anova_svm = Pipeline([('anova', anova_filter), ('svc', clf_kbest)])

###############################################################################
# Plot the cross-validation score as a function of percentile of features
score_means = list()
score_stds = list()
percentiles = (1, 3, 6, 10, 15, 20, 30, 40, 60, 80, 100)
kbest = (10, 20, 50, 100, 150, 200, 500, 1000, 5000 )

for k in kbest:
    print("kbest: " + str(k))
    anova_svm.set_params(anova__k=k, svc__C=.1)
    # Compute cross-validation score using all CPUs
    this_scores = cross_validation.cross_val_score(anova_svm, df_trans, sample_classification, n_jobs=1)
    score_means.append(this_scores.mean())
    score_stds.append(this_scores.std())

print("kbest:")
print(kbest)
print("score_means:")
print(score_means)
print("score_stds:")
print(score_stds)

plt.title(
    'Performance of the SVM-Anova varying the number of features selected')
plt.xlabel('k')
plt.ylabel('Prediction rate')

'''
for percentile in percentiles:
    print("percentile: " + str(percentile))
    clf.set_params(anova__percentile=percentile)
    # Compute cross-validation score using all CPUs
    this_scores = cross_validation.cross_val_score(clf, df_trans, sample_classification, n_jobs=1)
    score_means.append(this_scores.mean())
    score_stds.append(this_scores.std())

print("percentiles:")
print(percentiles)
print("score_means:")
print(score_means)
print("score_stds:")
print(score_stds)
plt.errorbar(percentiles, score_means, np.array(score_stds))

plt.title(
    'Performance of the SVM-Anova varying the percentile of features selected')
plt.xlabel('Percentile')
plt.ylabel('Prediction rate')
'''

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