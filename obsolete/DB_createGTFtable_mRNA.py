__author__ = 'KATRINA'



'''

Read in the GTF file and add all entries to the Gene table within the dropseq database

note: this was used to create the `Genes` table in initial dropseq database, but this table
has been removed from the dropseq2 database and a table containing only gene names exists

this script may be useful for reconstruction of gene information from gene names

'''

import MySQLdb
import sys
from string import Template
import json

cnx = MySQLdb.connect("localhost","mysql_user","balamuthia","dropseq")
cursor = cnx.cursor()

add_gene_template = Template("INSERT INTO Gene"
              "(chr, feature, start, end, strand, official_geneID, Name, transcriptID) "
              "VALUES ('$chr', '$feature', '$start', '$end', '$strand', '$official_geneID', '$Name', '$transcript_id')")

gtf_file = open(sys.argv[1],'r')
gtf_dictionary = {}

while True:
    x = gtf_file.readline()
    x = x.rstrip()
    if not x: break
    if x[0] != '#':
        info = x.split('\t')
        gene_info = info[8].split(';')
        input = add_gene_template.substitute(chr=info[0],feature=info[2], start = info[3],end = info[4], strand = info[6],official_geneID=gene_info[0].split()[1].replace("\"",""),Name=gene_info[1].split()[1].replace("\"","").strip("\""),transcript_id=gene_info[2].split()[1].replace("\"",""))
        cursor.execute(input)
        gene_id = cursor.lastrowid
        gtf_dictionary[gene_info[1].split()[1].replace("\"","").strip("\"")] = gene_id


cnx.commit()
cursor.close()
cnx.close()

json.dump(gtf_dictionary, open("/data/katrina/Scripts/gene_dictionary.txt",'w'))