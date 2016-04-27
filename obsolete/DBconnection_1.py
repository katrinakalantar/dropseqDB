__author__ = 'KATRINA'

'''

initial attempt at a direct connection to DROPSEQ DB
note: this script is extremely slow!

DBconnection_2.py [input_directory *containing all the sample files to be loaded with that experiment]

#outdated

'''


import MySQLdb
import os
import glob
from pandas import DataFrame
import sys
import random
import string
from string import Template

cnx = MySQLdb.connect("localhost","mysql_user","balamuthia","dropseq")
cursor = cnx.cursor()

add_experiment = ("INSERT INTO Experiment "
               "(Name, Description) "
               "VALUES (%s, %s)")

get_experimentID_template = Template("SELECT id FROM Experiment WHERE Name='$experiment_name'")

add_sample = ("INSERT INTO Sample"
              "(ExperimentID, SampleType, TissueType, Name) "
              "VALUES (%(exp_id)s, %(samp_type)s, %(tissue_type)s, %(name)s)")

add_sample_template = Template("INSERT INTO Sample"
              "(ExperimentID, SampleType, TissueType, Name) "
              "VALUES ('$exp_id', '$samp_type', '$tissue_type', '$name')")

get_sampleID_template = Template("SELECT id FROM Sample WHERE Name='$sample_name'")

'''
add_gene = ("INSERT INTO Gene"
              "(Accession, Name, Description) "
              "VALUES (%(gene_acc)s, %(gene_name)s, %(gene_desc)s")
add_gene_template = Template("INSERT INTO Gene"
              "(Accession, Name, Description) "
              "VALUES ('$gene_acc', '$gene_name', '$gene_desc')")
'''
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

add_data = ("INSERT INTO Data"
              "(SampleID, GeneID, CellID, ReadCounts) "
              "VALUES (%(sample_id)s, %(gene_id)s, %(cell_id)s, %(rc)s")

add_data_template = Template("INSERT INTO Data"
              "(SampleID, GeneID, CellID, ReadCounts) "
              "VALUES ('$sample_id', '$gene_id', '$cell_id', '$rc')")

def getGeneAccession(g):
    #CONNECT TO NCBI TO GET THIS INFO
    #TEMPORARY - until I actually figure out how to get the accession numbers from online DB
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(7))

def getGeneDescription(g):
    #CONNECT TO GENECARDS TO GET THIS INFO
    a = 1

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

    #if the sampleID is not in there, then add it
    if(sample_id == 0):

        df = DataFrame.from_csv(f, sep="\t")

        #ADD ALL GENES FROM FILE
        #genes2database(df.index)

        #ADD SAMPLE, CELLS, DATA
        print(sample_name)
        samp_type = "RNA"
        tissue_type = "CSF"
        input = add_sample_template.substitute(exp_id=exp_id,samp_type=samp_type, tissue_type=tissue_type,name=sample_name)
        cursor.execute(input)
        sample_id = cursor.lastrowid
        #print(sample_id)
    #else, continue on to add barcodes to exisiting sapmleID

        for barcode in df.columns.values: #loop through all barcodes
            #print(sample_id)
            #print(barcode)
            input = add_cell_without_annotation_template.substitute(sample_id=sample_id,barcode=barcode)
            cursor.execute(input)
            cell_id = cursor.lastrowid
            for i in range(len(df[barcode])): #loop through all genes within cell column, populate data
                gene_name = df.index[i]
                #print(gene_id)
                #print(gene.index)
                #print(df.at[i,barcode])
                #print(df.at[gene.index,barcode])
                read_count = df.at[gene_name,barcode]#GET READ COUNT
                #print(read_count)

                gene_query = get_geneID_template.substitute(gene_name=gene_name)
                gene_id = cursor.execute(gene_query)
                if(gene_id == 0): #gene does not yet exist in Gene DB, add it.
                    print("ERROR: gene " + str(gene_name) + " does not exist in the Gene table for HG38")
                else: #gene already exists, go ahead and add it
                    input = add_data_template.substitute(sample_id=sample_id, gene_id=gene_id, cell_id=cell_id, rc=read_count)
                    cursor.execute(input)

        #commit everything from that sample file
        cnx.commit()
        print("completed addition of sample " + sample_name + " to database")

#close cursor and connection after all samples have been added to the table
cursor.close()
cnx.close()

