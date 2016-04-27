__author__ = 'KATRINA'

'''
upload_data1.py

generate a .sql file which can then be executed within MySQL to upload the data

input:
 1. the input directory (containing counts files)
 2. the starting sampleID (will increment by 1 for every remaining sample file in the input directory)
 3. experimentID
 4. tissue type
 5. sample type (if not RNA)
a semi-manual / semi-automated method for uploading data - creates .sql files to upload cell information for a sample

'''

#import MySQLdb
import pymysql
import os
import glob
from pandas import DataFrame
import sys
import numpy as np

#normalize all cells in the data matrix to RPM per cell
def normalizePerIndividualCell(array_expn_values):
    normalized_data = []
    total_transcripts_per_cell = sum(array_expn_values)
    per_million_scaling_factor = float(total_transcripts_per_cell)/float(1000000)
    newArray = array_expn_values/per_million_scaling_factor
    return newArray

#set all input parameters based on the command line arguments
input_directory = sys.argv[1]
print("input_directory: " + input_directory)
os.chdir(input_directory)
tissue_type = sys.argv[4]
exp_id = sys.argv[3]
sample_id = int(sys.argv[2])
samp_type = "RNA"
if len(sys.argv)>5:
    samp_type = sys.argv[5]

#BEGIN ADDING SAMPLES, CELLS, DATA
filenames = glob.glob('./*.txt')
for f in filenames:
    print("filename: " + f)
    output_file = open(os.path.join(input_directory,f+'output.sql'),'w')
    sample_name = os.path.split(f)[-1].split('.')[0]
    df = DataFrame.from_csv(f, sep="\t",index_col=0)

    input = ("INSERT INTO Sample (ExperimentID, SampleType, TissueType, Name) VALUES ("+str(exp_id)+',\''+samp_type+'\', \''+tissue_type+'\', \'' + sample_name+'\')')
    print(input)

    if True:
        print("nothing found in gene_set")
        print("INSERT INTO GeneSet (GeneList) VALUES (''"+'\',\''.join(str(x) for x in df.index.values)+"'')")
        print(input)

    for barcode in df.columns.values: #loop through all barcodes (cells)
        raw_expn_data = np.array(df[barcode])
        orig_expn = sum(raw_expn_data)
        normalized_expn_data = normalizePerIndividualCell(raw_expn_data)
        gene_set_id=1
        input2 = "INSERT INTO Cell (SampleID, GeneSetID, Barcode, OriginalExpression, NormalizedData) VALUES ("+str(sample_id)+\
                                                        ", " + str(gene_set_id) + ", '"+barcode+"', " + str(orig_expn) + ", '" + \
                 (', '.join(str(x) for x in normalized_expn_data))+"');"
        output_file.write(input2+"\n")

    sample_id += 1 #increment sampleID for the next file iteration
