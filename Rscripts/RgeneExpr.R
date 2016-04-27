#source("https://bioconductor.org/biocLite.R")
#biocLite()
#source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")
library(DESeq); library(statmod); library(pcaMethods); library(fastICA);
require(DESeq)

#d1 <- t(read.table("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\Lupus_data\\P-19_S9_L007_R1_001_matrix.txt", header=TRUE, row.names=1, sep="\t"))
#colnames(d1) <- paste("P", colnames(d1), sep = "_")
#d2 <- t(read.table("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\Lupus_data\\PF-25_S13_L007_R1_001_matrix.txt", header=TRUE, row.names=1, sep="\t"))
#colnames(d2) <- paste("PF", colnames(d2), sep = "_")

d1 <- t(read.table("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\DerisiLab\\Data\\Lupus_data\\analysis_testipython\\smallerDF.txt", header=TRUE, row.names=1, sep="\t"))


#data1 <- merge(d1, d2, by=0, all=TRUE) #merge two dataframes by row names
data1 <- d1
means <- rowMeans(data1)
vars <- apply(data1, 1, var)
cv2 <- vars/means^2
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9)
smoothScatter(log(means),log(cv2))

biocLite("statmod")
library(statmod);
require(statmod)
minMeanForFit <- unname( quantile( means[ which( cv2 > .3 ) ], .95 ) )
useForFit <- means >= minMeanForFit # & spikeins
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
fit$coefficients


# repeat previous plot
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2));
xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
vfit <- a1/xg + a0
# add fit line
lines( log(xg), log(vfit), col="black", lwd=3 )
df <- ncol(data1) - 1
# add confidence interval
lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")


afit <- a1/means+a0
varFitRatio <- vars/(afit*means^2)
varorder <- order(varFitRatio,decreasing=T)
oed <- data1[varorder,]
# save for the next exercise
save(oed,file="oed.RData")

# repeat previous plot
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2)); lines( log(xg), log(vfit), col="black", lwd=3 ); lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black"); lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black");
# add top 100 genes
points(log(means[varorder[1:100]]),log(cv2[varorder[1:100]]),col=2)


pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
adj.pval <- p.adjust(pval,"fdr")
sigVariedGenes <- adj.pval<1e-3;
table(sigVariedGenes)



m <- oed[1:50,]
heatmap(m/apply(m,1,max),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",ColSideColors=ifelse(grepl("ES",colnames(m)),"red","blue"))

winsorize <- function (x, fraction=0.05) {
  if(length(fraction) != 1 || fraction < 0 ||
     fraction > 0.5) {
    stop("bad value for 'fraction'")
  }
  lim <- quantile(x, probs=c(fraction, 1-fraction))
  x[ x < lim[1] ] <- lim[1]
  x[ x > lim[2] ] <- lim[2]
  x
}

# winsorize to remove 2 most extreme cells (from each side)
wed <- t(apply(data1, 1, winsorize, fraction=2/ncol(data1)))

# now let's recalculate the most variable genes with the winsorized matrix (wed)
means = rowMeans(wed); vars = apply(wed,1,var); cv2 <- vars/means^2
useForFit <- means >= unname( quantile( means[ which( cv2 > .3 ) ], .95 ) ) 
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
afit <- fit$coef["a1tilde"]/means+fit$coef["a0"]
vfit <- fit$coef["a1tilde"]/xg+fit$coef["a0"]
varFitRatio <- vars/(afit*means^2)
varorder <- order(varFitRatio,decreasing=T)
oed <- wed[varorder,]
# save for the next exercise
save(oed,file="oed_win.RData")

xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2)); lines( log(xg), log(vfit), col="black", lwd=3 ); lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black"); lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black");
# add top 100 genes
points(log(means[varorder[1:100]]),log(cv2[varorder[1:100]]),col=2)


m <- oed[1:50,]
heatmap(m/apply(m,1,max),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",ColSideColors=ifelse(grepl("ES",colnames(m)),"red","blue"))

plot(hclust(as.dist(1-cor(data1))))
