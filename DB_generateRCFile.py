__author__ = 'KATRINA'

'''
March 18, 2016

python DB_generateRCFile.py "SELECT * FROM Data WHERE (SampleID=3 OR SampleID=5) AND (GeneID=4 OR GeneID=10)" outputFile [-RPM]

input: query to DB
output: resulting counts file

TODO:
add a dictionary lookup from id to name, so files contain geneNames instead of IDs
'''

import MySQLdb
import sys
import numpy as np
import pandas as pd
import json


#convert input (Raw) dataframe into RPM normalized read count file
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


print(sys.argv[1])

cnx = MySQLdb.connect("localhost","mysql_user","balamuthia","dropseq")
cursor = cnx.cursor()

query_partA = sys.argv[1].split('*')[0]
query_partB = sys.argv[1].split('*')[1]
print(query_partA)
print(query_partB)


get_unique_geneIDs = query_partA + "GeneID" + query_partB #create query to get all GeneIDs for the set
cursor.execute(get_unique_geneIDs)
all_geneIDs = cursor.fetchall()
x = []
for i in all_geneIDs:
    x.append(i[0])
set_geneIDs = np.array(list(set(x)))

get_unique_cellIDs= query_partA + "CellID" + query_partB #create query to get all CellIDs for the set
cursor.execute(get_unique_cellIDs)
all_cellIDs = cursor.fetchall()
x = []
for i in all_cellIDs:
    x.append(i[0])
set_cellIDs = np.array(list(set(x)))

a = np.zeros(shape=(len(set_geneIDs),len(set_cellIDs)))
df = pd.DataFrame(a,columns=set_cellIDs, index=set_geneIDs) #create dataframe of 0's with col and row IDs


first_pass = True
sample_classification = []
get_data_query = sys.argv[1]   #run the input data table query
cursor.execute(get_data_query)
row = cursor.fetchone()
#iterate through rows of the query and add counts data to dataframe in corresponding position
#also, generate a sample classification based on sampleID
while row is not None:
    if first_pass:
        gid=row[2] #assuming this is the GeneID
        sample_classification.append(row[1])
        first_pass = False
    else:
        if row[2] == gid:
            sample_classification.append(row[1])
    df[row[3]][row[2]] = row[4]
    row = cursor.fetchone()

classification = {}
classification["SampleID"] = sample_classification
json.dump(classification, open(sys.argv[2]+"_classification.txt",'w'))


#print(df)
cursor.close()
cnx.close()

if "-RPM" in sys.argv:
    df_norm = normalizeRPM(df)
    df_norm.to_csv(sys.argv[2]+"_RPM.tsv",sep="\t",header=True,index=True)
else:
    df.to_csv(sys.argv[2],sep="\t",header=True,index=True) #save dataframe to a file (input file name)
