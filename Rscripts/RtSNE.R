#d1 <- t(read.table("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\hannas_dropseq\\N0_S1_R1_matrix.txt", header=TRUE, row.names=1, sep="\t"))
#d2 <- t(read.table("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\hannas_dropseq\\N1_S2_R1_matrix.txt", header=TRUE, row.names=1, sep="\t"))
#d3 <- t(read.table("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\hannas_dropseq\\N2_S3_R1_matrix.txt", header=TRUE, row.names=1, sep="\t"))
#d4 <- t(read.table("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\hannas_dropseq\\N3_S5_R1_matrix.txt", header=TRUE, row.names=1, sep="\t"))
#d5 <- t(read.table("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\hannas_dropseq\\N5_S7_R1_matrix.txt", header=TRUE, row.names=1, sep="\t"))
library(Rtsne)
data <- t(read.table("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\Lupus_data\\analysis_testipython\\PF25only\\genesubsetDF.txt", header=TRUE, row.names=1, sep="\t"))
print(data)
dataT <- t(as.matrix(data))

CAGIdata <- t(read.table("C:\\cygwin64\\home\\KATRINA\\UCB\\CourseMaterials\\COMPBIO290\\final_project\\CMPBIO290-2016_project\\katrinas_preliminary_analysis\\fasta_alt_DF.txt", header=TRUE, row.names=1, sep="\t"))
CAGI_unique = unique(CAGIdata)
#dataT <- t(as.matrix(CAGIdata))
#dim(dataT)
rtsne_out <- Rtsne(CAGI_unique)
plot(rtsne_out$Y, main="BarnesHutSNE",pch = 4, cex = .8)

