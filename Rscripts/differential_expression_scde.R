
#source("https://bioconductor.org/biocLite.R")
#biocLite("scde")
library("scde") #installed this via github: install_githup('hms-dbmi/scde',build_vignettes = FALSE)

# load dataset
dataset = "C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\Lupus_data\\analysis_testipython\\P19PF25_d3_top20\\smallerDF_scde.txt"
d1 <- as.matrix(read.table(dataset, header=TRUE, row.names=1, sep="\t"))

# factor determining cell types
sg <- factor(gsub("(X100|X200).*", "\\1", colnames(d1)), levels = c("X100","X200"))

# the group factor should be named accordingly
names(sg) <- colnames(d1)
table(sg)

# clean up the dataset
cd <- clean.counts(d1, min.lib.size=1000, min.reads = 1, min.detected = 1) #this didn't seem to work for d1 or m

# reg-executing this on cd so that sg and cd are same legth for the subsequent section
sg <- factor(gsub("(X100|X200).*", "\\1", colnames(cd)), levels = c("X100","X200"))


# EVALUATION NOT NEEDED - THIS TAKES A LONG TIME SUPPOSEDLY ON JUST 1 CORE...BE PREPARED
# calculate models
#counts<-apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)


###### WAITING FOR STEP TO COMPLETE.

# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- o.ifm$corr.a > 0
table(valid.cells)

o.ifm <- o.ifm[valid.cells, ]

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = TRUE)




# define two groups of cells
groups <- factor(gsub("(X100|X200).*", "\\1", rownames(o.ifm)), levels  =  c("X100", "X200"))
names(groups) <- row.names(o.ifm)
# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

######running at 3:39 pm

# top upregulated genes (tail would show top downregulated ones)
head(ediff[order(ediff$Z, decreasing  =  TRUE), ])

# write out a table with all the results, showing most significantly different genes (in both directions) on top
output_file = "C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\Lupus_data\\analysis_testipython\\P19PF25_d3_top20\\DEresults.txt"
write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = output_file, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)



#alternatively, we can run DE analysis on a single gene and visualize the results
scde.test.gene.expression.difference("IL7R", models = o.ifm, counts = cd, prior = o.prior)



# get failure probabilities on the expresison range
o.fail.curves <- scde.failure.probability(o.ifm, magnitudes = log((10^o.prior$x)-1))
par(mfrow = c(1,1), mar = c(3.5,3.5,0.5,0.5), mgp = c(2.0,0.65,0), cex = 1)
plot(c(), c(), xlim=range(o.prior$x), ylim=c(0,1), xlab="expression magnitude (log10)", ylab="drop-out probability")
invisible(apply(o.fail.curves[, grep("X100",colnames(o.fail.curves))], 2, function(y) lines(x = o.prior$x, y = y,col = "orange")))
invisible(apply(o.fail.curves[, grep("X200", colnames(o.fail.curves))], 2, function(y) lines(x = o.prior$x, y = y, col = "dodgerblue")))

# get failure probabilities on the expresison range
o.fail.curves <- scde.failure.probability(o.ifm, magnitudes = log((10^o.prior$x)-1))
# get self-fail probabilities (at a given observed count)
p.self.fail <- scde.failure.probability(models = o.ifm, counts = cd)



###PAGODA
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)



