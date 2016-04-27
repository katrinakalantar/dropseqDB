# DropSeqViewer 
Single Cell RNA-seq data visualization tool, database, and data manipulation scripts
work done in the Derisi Lab at UCSF (Spring 2016 Rotation)

#### Most Important Scripts:
##### upload_data1.py
##### d3pipeline.py
##### create_dropseqDB2.sql
##### dropseq_1.py
##### normalizeRPM.py

### CREATE DROPSEQ2 DATABASE

#### create_dropseqDB2
This file contains the .sql specification for the current implementation of the DROPSEQ2 database


### UPLOAD DATA FROM FILES TO DATABASE

####upload_data1.py
generate a .sql file which can then be executed within MySQL to upload the data
input:
 1. the input directory (containing counts files)
 2. the starting sampleID (will increment by 1 for every remaining sample file in the input directory)
 3. experimentID
 4. tissue type
 5. sample type (if not RNA)
a semi-manual / semi-automated method for uploading data - creates .sql files to upload cell information for a sample


### WORKING WITH D3 VISUALIZATION

#### d3pipeline.py
from a counts file, classification file, and config file (specifying the method to use
for gene selection, thresholds, and dimensionality reduction), create a d3data.tsv file which can then be
used with the DropSeqViewer. Also creates an _DFresult.txt file containing only the subset
of genes selected.
```
> python d3pipeline.py [input DF file] [input classification file] [config file]
```
example config file:
{
"gene_set_method":"kbest",
"gene_number":20,
"cell_expression_threshold":2000,
"dim_red_method":"ZIFA"
}

note:
- input dataframe is expected to have cells in the column orientation (and genes in rows)
- limited error checking on the config file, if you are missing a field in the config file
the script will likely error. TODO


### OTHER BASIC USE FILES
#### dropseq_1.py
Pipeline-style script containing many functions to remove low expression rows/columns (genes/cells)
and remove high variance data (genes/cells).
Saves the final dataframe as [input_file]_filtered.txt
```
> python dropseq_1.py [input_file]
```
note: this script has hard-coded threshold for removal of data. It is advised that you ovserve the distributions
within a particular sample using investigateGeneVariance.py, plotGeneExprDistribution.py etc prior to running
this script and modify the parameters accordingly.


####normalizeRPM.py
Normalize an input file (hard-coded) 
This script contains both the original method for normalization (full file RPM), which
is actually not correct.
AND the updated normalization method "normalizeRPM_new" which normalizes per cell (column)


###//investigative
goal: OBSERVE DISTRIBUTIONS WITHIN THE DATA:

####investigateGeneVariance.py
Plot the distribution of variance across genes and cells.
```
> python investigateGeneVariance.py [file of interest]
```

####plotGeneExprDistribution.py

plot the distributions of transcript counts across
1. genes
2. cells
3. specific gene across all cells
```
>python plotGeneExprDistribution.py [input_file]
```

####datasetMetrics_1.py
generate a data frame with statistics on expression for each gene
options:
1. sort on total expression summed over the whole sample (all cells)
2. sort on mean expression for only cells expressing that gene
```
> python datasetMetrics_1.py [input file]
```


####clustermap_1.py
Generate a hierarchical clustering heatmap of the input dataframe
```
> python clustermap_1.py [input_file]
```
note: in a sparse matrix this will take a long time to execute and provide
minimal useful information. Best to used this on a sub-selected matrix.


####heatmap.py
create a heatmap from input dataframe
note: I was using this specifcially with output of BackSPIN algorithm
```
> python heatmap.py [input_file]
```

###//DropseqViewer
goal: this directory contains all files required to load webserver, including a sample data file

app.js, index.html, and style.css are the primary files
however, they rely on FileSaver.js, require.js, script.js, and summary.html for full loading

the /analysis directory contains an example_lupusPvNID dataset with a d3data.tsv file to test loading


###//obsolete
note: these files are outdated or no longer useful in the direct analysis pipeline

####CSFTest_1.py
Playing around with CSF data -
loads all data into a dataframe and creates the classification vector
the combined data is saved in a CSF_combined_matrix.txt file and the classification vector
is saved in classification.json


####d3pipeline_2bin.py
same as the original d3pipeline_1.py, but this one transforms the dataframe into a binary frame prior to apply tsne...attempting to see if that will solve the variance issue
```
> python d3pipeline.py [input DF file] [input classification file] [config file]
```


####DB_createGTFtable.py
Read in the GTF file and add all entries to the Gene table within the dropseq database
note: this was used to create the `Genes` table in initial dropseq database, but this table
has been removed from the dropseq2 database and a table containing only gene names exists

this script may be useful for reconstruction of gene information from gene names


####DB_createGTFtable_mRNA.py
Read in the GTF file and add all entries to the Gene table within the dropseq database

note: this was used to create the `Genes` table in initial dropseq database, but this table
has been removed from the dropseq2 database and a table containing only gene names exists

this script may be useful for reconstruction of gene information from gene names


####DB_generateRCFile.py
NOTE: this file was for reconstruction from initial dropseq DB to flat file, it is not useful with DROPSEQ2 DB
```
> python DB_generateRCFile.py "SELECT * FROM Data WHERE (SampleID=3 OR SampleID=5) AND (GeneID=4 OR GeneID=10)" outputFile [-RPM]
```
input: query to DB
output: resulting counts file

TODO:
add a dictionary lookup from id to name, so files contain geneNames instead of IDs


####DBconnection_1.py
initial attempt at a direct connection to DROPSEQ DB
note: this script is extremely slow!
```
> python DBconnection_2.py [input_directory *containing all the sample files to be loaded with that experiment]
```


####DBconnection_2.py
this file differs from DBconnection_1.py because it creates a new file format from the table and then loads that into the mysql database
```
> python DBconnection_2.py [input_directory *containing all the sample files to be loaded with that experiment]
```
initial attempt at a direct connection to DROPSEQ DB
note: this script is [ALSO] extremely slow!


####DBconnection_3.py
```
> python DBconnection_3.py [path to experiment dir containing read count files]
```
add the experiment, samples, cell, and individual data count values to the database
updates will only be committed if the whole sample file is updated without interruption

note: only use with initial DROPSEQ DB;
this version WAS functional, but large datasets took forever to load


####DBconnection_4.py
Another outdated version of database upload script - actually unsure what this one was for
or how it was different from the previous 3 


####GenerateD3File_1.py
```
> python GenerateD3File_1.py [input dataframe] [input classification .json file] [-o=output_file_name] [-s=[GENE1,GENE2,GENE3]]
```
note: this was an early version preceding d3pipeline.py


####GenerateD3File_2binary.py
```
> python GenerateD3File_2.py [input dataframe] [input classification .json file] [-o=output_file_name] [-s=[GENE1,GENE2,GENE3]]
same as the original GenerateD3File_1.py, but this one transforms the dataframe into a binary frame prior to apply tsne...attempting to see if that will solve the variance issue
```
note: this was an early version preceeding d3pipeline_2bin.py


####ML_1.py
investigating Machine Learning (ML) methods for determining the most informative genes for separating two datasets
```
> python ML_1.py inputDF classificationFile [-norm]
```

####LupusTest_local.py, LupusTest_1.py
This script was used to investigate the data associated with lupus v healthy. It was not intended to 
perform a particular task or return a particular result, but rather was continually modified to
try to ways of manipulating the data


####topGenes_1.py
generate a data frame with statistics on expression for each gene
options:
1. sort on total expression summed over the whole sample (all cells)
2. sort on mean expression for only cells expressing that gene

note: this script was used to determine the most highly expressed genes from CSF data


####tSNE_1.py
first attempt at implementing TSNE
http://scikit-learn.org/stable/auto_examples/manifold/plot_lle_digits.html

####tSNE_2.py
Create Combined sample TSNE plot

1. Loop through input folder
2. Gather all .txt files (all .txt files in folder should correspond to counts tables, containing the sample genes)
3. Create a combined dataframe for all samples in the folder
    - keep track of where the data frame entries came from
4. Clean the data frame with filters on rows and columns
5. Run TSNE algorithm
6. Plot the TSNE with color corresponding to which sample they came from

note: this was prohibitively slow on my laptop - I saved the combined file and ran TSNE in R


####ZIFA_test1.py
Run ZIFA algorithm on input file, print the result; this is intended to test the method. 
Implementation of the result was done in ipython notebook.
```
> python ZIFA_test1.py [input_file]
```