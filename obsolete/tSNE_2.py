__author__ = 'KATRINA'

'''

Create Combined sample TSNE plot

1. Loop through input folder
2. Gather all .txt files (all .txt files in folder should correspond to counts tables, containing the sample genes)
3. Create a combined dataframe for all samples in the folder
    - keep track of where the data frame entries came from
4. Clean the data frame with filters on rows and columns
5. Run TSNE algorithm
6. Plot the TSNE with color corresponding to which sample they came from

note: this was prohibitively slow on my laptop - I saved the combined file and ran TSNE in R

'''

from pandas import DataFrame
import os
import glob

def removeRowsLessThanTen(df):
    df_no_under_ten = df.loc[(df.sum(axis=1) >= 10)] #remove all rows that have sum < 10
    return df_no_under_ten

def removeCols(df):
    df_no_col_under_twothousand = df.loc[:,(df.sum(axis=0) >= 2000)]
    return df_no_col_under_twothousand


#set input dir, move there
input_directory = "C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\test_data\\CSF"
os.chdir(input_directory)
filenames = glob.glob('./*.txt')

df_merged = DataFrame
labels_merged = []
count = 0

#loop through files in input dir
for f in filenames:
    count += 1
    print(f)
    df = removeCols(DataFrame.from_csv(f, sep="\t")) #remove cols with total expn < 2000 counts initially
    a = [count for i in range(len(df.columns.values))] #create labels array to classify this set of cells in final df

    try: #try to merge the data frame to the existing df_merged data frame
        dfm = df_merged
        df_merged = dfm.merge(df, left_index=True, right_index=True, how = 'right')
    except TypeError: #if error because df_merged is uninitialized, set df_merged to df
        df_merged = df

    labels_temp = labels_merged #merge the labels array from this iteration to the total labels array
    labels_merged = labels_temp + a

    print("original file #cells: " + str(len(df.columns.values)))   #print statistics on original file
    print("original file #genes: " + str(len(df.index)))
    print("merged file #cells: " + str(len(df_merged.columns.values)))   #print statistics on merged file
    print("merged file #genes: " + str(len(df_merged.index)))

df_clean = removeRowsLessThanTen(df_merged) #remove all gene rows with expression less than threshold
print(labels_merged)

print("filtered file #cells: " + str(len(df_clean.columns.values)))   #print statistics on what was filtered out
print("filtered file #genes: " + str(len(df_clean.index)))

'''

#THIS PART IS TOO SLOW TO RUN ON MY PERSONAL COMPUTER

#Initialize variables to run TSNE
X = df_clean
n_samples, n_features = X.shape[0],X.shape[1]
n_neighbors = 30

print("Computing t-SNE embedding")
tsne = manifold.TSNE(n_components=2, init='pca', random_state=0)
X_tsne = tsne.fit_transform(X)

vis_x = X_tsne[:, 0]
vis_y = X_tsne[:, 1]
y_data = labels_merged

plt.scatter(vis_x, vis_y,c=y_data, cmap=plt.cm.get_cmap("jet", count))
plt.colorbar(ticks=range(count))
plt.clim(-0.5, 9.5)
plt.show()

'''

#INSTEAD, ON MY LAPTOP...TRANSPOSE AND SAVE TO LOAD INTO RtSNE
(df_clean.transpose()).to_csv("full_tSNE_data.tsv",sep="\t",header=True, index=True)
