__author__ = 'KATRINA'

'''
DBconnection_3.py [path to experiment dir containing read count files]

add the experiment, samples, cell, and individual data count values to the database
updates will only be committed if the whole sample file is updated without interruption

note: only use with initial DROPSEQ DB;
this version WAS functional, but large datasets took forever to load

#outdated

'''

import MySQLdb
import os
import glob
from pandas import DataFrame
import sys
from string import Template


cnx = MySQLdb.connect("localhost","mysql_user","balamuthia","dropseq")
cursor = cnx.cursor()

add_experiment = ("INSERT INTO Experiment "
               "(Name, Description) "
               "VALUES (%s, %s)")

get_experimentID_template = Template("SELECT id FROM Experiment WHERE Name='$experiment_name'")

add_sample_template = Template("INSERT INTO Sample"
              "(ExperimentID, SampleType, TissueType, Name) "
              "VALUES ('$exp_id', '$samp_type', '$tissue_type', '$name')")

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

add_data_template = Template("INSERT IGNORE INTO Data"
              "(SampleID, GeneID, CellID, ReadCounts) "
              "VALUES ('$sample_id', '$gene_id', '$cell_id', '$rc')")

input_directory = sys.argv[1]
print(input_directory)
os.chdir(input_directory)

experiment = os.path.split(input_directory)[-1]
experiment_query = get_experimentID_template.substitute(experiment_name=experiment)
exp_id = cursor.execute(experiment_query)
print(experiment)

if exp_id == 0: #create the new experiment in the database
    data_experiment = (str(experiment), 'first CSF experiment')
    cursor.execute(add_experiment, data_experiment)
    exp_id = cursor.lastrowid

print(exp_id)

#BEGIN ADDING SAMPLES, CELLS, DATA
filenames = glob.glob('./*.txt')
for f in filenames:
    print(f)
    sample_name = os.path.split(f)[-1].split('.')[0]
    sample_query = get_sampleID_template.substitute(sample_name=sample_name)
    sample_id = cursor.execute(sample_query)

    #if the sampleID is not in Samples table already, add it
    if(sample_id == 0):

        df = DataFrame.from_csv(f, sep="\t")

        #ADD SAMPLE, CELLS, DATA
        print(sample_name)
        samp_type = "RNA"
        tissue_type = "CSF"
        input = add_sample_template.substitute(exp_id=exp_id,samp_type=samp_type, tissue_type=tissue_type,name=sample_name)
        cursor.execute(input)
        sample_id = cursor.lastrowid

        for barcode in df.columns.values: #loop through all barcodes (cells)
            input = add_cell_without_annotation_template.substitute(sample_id=sample_id,barcode=barcode) #add cell to Cell table
            cursor.execute(input)
            cell_id = cursor.lastrowid
            for i in range(len(df[barcode])):        #loop through all genes within cell column, populate data
                gene_name = df.index[i]
                try:
                    read_count = df.at[gene_name,barcode]       #GET READ COUNT
                    gene_query = get_geneID_template.substitute(gene_name=gene_name)
                    cursor.execute(gene_query)
                    gene_id = cursor.fetchone()[0]
                    #print(gene_id)
                    input = add_data_template.substitute(sample_id=sample_id, gene_id=gene_id, cell_id=cell_id, rc=read_count)
                    cursor.execute(input)
                except:
                    print("Failed to add gene: " + str(gene_name))
            print('added cell ' + barcode)

        cnx.commit()    #commit all table updates from this sample file
        print("COMPLETED addition of sample " + sample_name + " to database")

#close cursor and connection after all samples have been added to the table
cursor.close()
cnx.close()

