__author__ = 'KATRINA'

'''

Annother outdated version of database upload script - actually unsure what this one was for
or how it was different from the previous 3

#outdated

'''

import MySQLdb
import os
import glob
from pandas import DataFrame
import pandas as pd
import sys
from string import Template
import json
import numpy as np

cnx = MySQLdb.connect("localhost","mysql_user","balamuthia","dropseq", local_infile=1)
cursor = cnx.cursor()

add_experiment = ("INSERT INTO Experiment "
               "(Name, Description) "
               "VALUES (%s, %s)")
get_experimentID_template = Template("SELECT id FROM Experiment WHERE Name='$experiment_name'")
add_gene_template_light = Template("INSERT INTO Gene"
              "(Name) "
              "VALUES ('$Name')")
add_sample = ("INSERT INTO Sample"
              "(ExperimentID, SampleType, TissueType, Name) "
              "VALUES (%(exp_id)s, %(samp_type)s, %(tissue_type)s, %(name)s)")
add_sample_template = Template("INSERT INTO Sample"
              "(ExperimentID, SampleType, TissueType, Name, TotalReads) "
              "VALUES ('$exp_id', '$samp_type', '$tissue_type', '$name', '$total_reads')")
get_sampleID_template = Template("SELECT id FROM Sample WHERE Name='$sample_name'")
get_geneID_template = Template("SELECT id FROM Gene WHERE Name='$gene_name'")
add_annotation_code = ("INSERT INTO AnnotationCodes"
              "(Annotation) "
              "VALUES (%(annot_code)s")
add_cell = ("INSERT INTO Cell"
              "(SampleID, AnnotationID, Barcode) "
              "VALUES (%(sample_id)s, %(annot_id)s, %(barcode)s")
add_cell_without_annotation_template = Template("INSERT INTO Cell"
              "(SampleID, Barcode) "
              "VALUES ('$sample_id', '$barcode')")
get_cellID_template = Template("SELECT id FROM Cell WHERE Barcode='$barcode' and SampleID='$sample_id'")
add_data = ("INSERT INTO Data"
              "(SampleID, GeneID, CellID, ReadCounts) "
              "VALUES (%(sample_id)s, %(gene_id)s, %(cell_id)s, %(rc)s")
get_data_template = Template("SELECT id FROM Data WHERE CellID='$cell_id' and GeneID='$gene_id'")
add_data_template = Template("INSERT IGNORE INTO Data"
              "(SampleID, GeneID, CellID, ReadCounts) "
              "VALUES ('$sample_id', '$gene_id', '$cell_id', '$rc')")

'''
OUTLINE
add experiment
query to get experiment ID
for every gene in the input file:
    query gene table to see if gene exists
    if gene does not exist in gene table
        add gene
add sample
query to get sampleID
for every cell in the input file:
    add cell to cell table
    query to get cellID
    for every gene in the input file:
        query to get geneID
        add data (sampleID, cellID, geneID)
'''

input_directory = sys.argv[1]
print(input_directory)
os.chdir(input_directory)

#load the gene dictionary from file
# ( file created by DB_createGTFtable_mRNA.py, run before this script to populate all
#   "expected" genes into the database )
gene_dict = json.load(open("/data/katrina/Scripts/gene_dictionary.txt"))


experiment = os.path.split(input_directory)[-1] # input directory name = ExperimentName
experiment_query = get_experimentID_template.substitute(experiment_name=experiment) # check if experiment exists in DB
cursor.execute(experiment_query)
x = cursor.fetchall()
if x == (): # nothing returned by query, no experiment exists by this name
    exp_id = 0
else:
    exp_id = x[0]

if exp_id == 0: # if experiment does not already exist, create the new experiment in the database
    data_experiment = (str(experiment), 'first CSF experiment')
    cursor.execute(add_experiment, data_experiment)
    exp_id = cursor.lastrowid

#BEGIN ADDING SAMPLES, CELLS, DATA
filenames = glob.glob('./*.txt')
config = json.load(open('experiment.config','r')) # hard-coded config file, required to run
for f in filenames: # loop over all .txt files in the directory (assumes that only counts tables are in directory)
    print(f)
    sample_name = os.path.split(f)[-1].split('.')[0]

    sample_query = get_sampleID_template.substitute(sample_name=sample_name)
    cursor.execute(sample_query) # check if experiment exists in DB
    x = cursor.fetchall()
    if x == (): # nothing returned, sample did not previously exist in DB
        sample_id = 0
    else:
        sample_id = x[0][0]

    output_file=open(sample_name+"out.tsv",'w') # create a file containing Data values to load into DB manually

    # calculated the total reads within a file
    total_transcripts_per_cell = []
    df = DataFrame.from_csv(f, sep="\t", index_col=0) # read in the counts table
    for i in df.columns.values:
        sum = 0
        for genes in df[i]:
            sum += genes
        total_transcripts_per_cell.append(sum)
    total_sum = np.sum(total_transcripts_per_cell)

    if(sample_id == 0): # if the sampleID is not in the DB, then add it
        samp_type = config["sample_type"]    # get the sample type from the config file
        tissue_type = config["tissue_type"]     # get the tissue type from the config file
        try:
            x = exp_id[0]
        except:
            x=exp_id
        exp_id = x
        input = add_sample_template.substitute(exp_id=exp_id,samp_type=samp_type,tissue_type=tissue_type,name=sample_name, total_reads=total_sum)
        cursor.execute(input)
        sample_id = cursor.lastrowid


    count = 0
    output = [] #create the array which will be converted to mysql df upload

    for i in df.columns.values: # loop through all cells (columns in DF)
        count +=1
        print("BEGIN NEW CELL: " + str(count) + " of " + str(len(df.columns.values)))
        cell_query = get_cellID_template.substitute(sample_id=sample_id, barcode=i)
        cursor.execute(cell_query) # check if the cell id exists in the Cell table
        x = cursor.fetchall()
        cell_id=0
        if x == (): # cell does not exist, create it
            input = add_cell_without_annotation_template.substitute(sample_id=sample_id,barcode=i)#[0],barcode=i)
            cursor.execute(input)
            cell_id = cursor.lastrowid
        else: #cell exists, use the existing cellID
            cell_id = x[0]

        for j in list(df.index): # loop through all geneIDs (rows in DF)
            try:
                gene_id = gene_dict[j] # try to get the gene_id from the gene_dict
            except: # the geneID did not exist in in the Gene table, add it and update gene_id
                if str(j) not in ["","nan"]: # make sure the gene_id is not "nan"
                    print("ERROR: gene " + str(j) + " does not exist in the Gene table for HG38")
                    add_gene_light_query = add_gene_template_light.substitute(Name=str(j))
                    cursor.execute(add_gene_light_query)
                    gene_id = cursor.lastrowid
                    gene_dict[j] = gene_id

            try:
                data_query = get_data_template.substitute(cell_id=cell_id[0],gene_id=gene_id)
            except:
                data_query = get_data_template.substitute(cell_id=cell_id,gene_id=gene_id)
            cursor.execute(data_query) # check to see if Data entry exists in the DB
            d = cursor.fetchall()
            if d == (): #data does not exist, add it to the file of data to load manually.
                if str(j) not in ["","nan"]: # make sure the gene_id is not "nan"
                    try:
                        arr_line = ','.join([str(sample_id),str(gene_id),str(cell_id[0]),str(df[i][j])])
                        #input = add_data_template.substitute(sample_id=sample_id, gene_id=gene_id, cell_id=cell_id[0], rc=df[i][j])
                        #cursor.execute(input)
                    except:
                        arr_line = ','.join([str(sample_id),str(gene_id),str(cell_id),str(df[i][j])])
                        #input = add_data_template.substitute(sample_id=sample_id, gene_id=gene_id, cell_id=cell_id, rc=df[i][j])
                        #cursor.execute(input)
                    #output_file.writelines(arr_line + "\n")
                    output.append(arr_line)

        cnx.commit # commit each cell once it is complete

    uploadDF= pd.DataFrame(output, columns=['Sample_ID', 'GeneID','CellID','ReadCounts'])
    uploadDF.to_sql('Data',cnx,if_exists='append',index=False,chunksize=1000)



    #cnx.commit()
    #print(sample_id)
    ###TODO: Having difficulty with this - I am just creating the files (in the folder containing the counts tables) now and will load them manually after completion.
    #load_data_template = Template("LOAD DATA LOCAL INFILE '$out_file' INTO TABLE Data FIELDS TERMINATED BY ',' (SampleID,GeneID,CellID,ReadCounts);")# SET id = NULL;"# SET id = NULL;"
    #load_data_query = load_data_template.substitute(out_file=os.path.join(sys.argv[1],'out.tsv')) #'/data/katrina/Data/CSF_Experiment1/out.tsv')#
    #cursor.execute(load_data_query)
    #print('loaded data')
    output_file.close()


    #commit everything from that sample file
    cnx.commit()
    print("completed addition of sample " + sample_name + " to database (sans data)")

#close cursor and connection after all samples have been added to the table
cursor.close()
cnx.close()
json.dump(gene_dict, open("/data/katrina/Scripts/gene_dictionary.txt",'w'))
