__author__ = 'KATRINA'

'''
DBconnection_DS2.py [path to experiment dir containing read count files]

add the experiment, samples, cell, and individual data count values to the database
updates will only be committed if the whole sample file is updated without interruption

#NOTE: THIS IS INTENDED FOR USE WITH DROPSEQ2 DB, however I was unable to get the
connection (cnx) working with any of the imported python libraries despite previously
being able to connect to DROPSEQ DB with pymsql.
'''

#import MySQLdb
import pymysql
import mysql.connector
import os
#from sqlalchemy.dialects.mysql import pymysql
import glob
from pandas import DataFrame
import sys
import numpy as np
from sqlalchemy import create_engine
from string import Template

#cnx = mysqldb.connect("localhost","mysql_user","balamuthia","dropseq2")
#cnx = pymysql.connect(user='mysql_user', passwd = 'balamuthia', database='dropseq2',host='localhost', port=3306)#, unix_socket="/var/run/mysqld/mysqld.sock")#'169.230.81.133') #host='127.0.0.1',
#engine = create_engine("mysql+pymysql://mysql_user:balamuthia@localhost/db?host=localhost?port=3306")
#cnx = engine.connect()
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

get_geneSet_query = ("SELECT * FROM GeneSet")
add_geneSet_template = Template("INSERT INTO GeneSet"
                                "(GeneList)"
                                "VALUES (%(gene_list)s")

add_annotation_code = ("INSERT INTO AnnotationCodes"
              "(Annotation) "
              "VALUES (%(annot_code)s")

add_cell = ("INSERT INTO Cell"
              "(SampleID, AnnotationID, Barcode, OriginalExpression, ExpressionData) "
              "VALUES (%(sample_id)s, %(annot_id)s, %(barcode)s, %(orig_expn)s,%(expn_data)s")

add_cell_without_annotation_template = Template("INSERT INTO Cell"
              "(SampleID, GeneSetID, Barcode, OriginalExpression, NormalizedData) "
              "VALUES (%(sample_id)s,%(gene_set_id)s, %(barcode)s, %(orig_expn)s,%(norm_data)s")

add_data_template = Template("INSERT IGNORE INTO Data"
              "(SampleID, GeneID, CellID, ReadCounts) "
              "VALUES ('$sample_id', '$gene_id', '$cell_id', '$rc')")


def normalizePerIndividualCell(array_expn_values):
    print("normalizing input data per cell")
    normalized_data = []
    total_transcripts_per_cell = sum(array_expn_values)
    per_million_scaling_factor = float(total_transcripts_per_cell)/float(1000000)
    newArray = array_expn_values/per_million_scaling_factor
    print(sum(newArray))
    return newArray


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

        df = DataFrame.from_csv(f, sep="\t",index_col=0)

        #ADD SAMPLE, CELLS, DATA
        print(sample_name)
        samp_type = "RNA"
        tissue_type = "CSF"
        input = add_sample_template.substitute(exp_id=exp_id,samp_type=samp_type, tissue_type=tissue_type,name=sample_name)
        cursor.execute(input)
        sample_id = cursor.lastrowid

        cursor.execute(get_geneSet_query)
        gene_set = cursor.fetchall()
        print(gene_set)
        if len(gene_set)<1: #nothing in the gene set, add the current gene set.
            print("nothing found in gene_set")
            input = add_geneSet_template.substitute(gene_list=str(df.index.values))
            cursor.execute(input)
            gene_set_id = cursor.lastrowid
        else:
            gene_set_id = 0
            for x in gene_set: #loop through gene sets
                print(x)
                #check if the gene set matches str(df.index.values)
                #else, add the gene set

        for barcode in df.columns.values: #loop through all barcodes (cells)

            raw_expn_data = np.array(df[barcode])

            orig_expn = sum(raw_expn_data)
            normalized_expn_data = normalizePerIndividualCell(raw_expn_data)

            input = add_cell_without_annotation_template.substitute(sample_id=sample_id, gene_set_id=gene_set_id, barcode=barcode,
                                                                    orig_expn=orig_expn, norm_data=normalized_expn_data) #add cell to Cell table
            cursor.execute(input)
            cell_id = cursor.lastrowid

            print('added cell ' + barcode)

        cnx.commit()    #commit all table updates from this sample file
        print("COMPLETED addition of sample " + sample_name + " to database")

#close cursor and connection after all samples have been added to the table
cursor.close()
cnx.close()

