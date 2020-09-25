# mark individual genes on histogram of Stan results

library(lme4)
library(arm)
library(rstan)
library(coda)
library(gtools)
library(ggplot2)
library(bayesplot)
library(reshape2)

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# Yanlin's working directory
setwd("D:/DataPlus2017/Data")

# load annotation, Stan results, and gene scores
load("myFA_bulkonly.Rdata")
#load("StanCfullResults.Rdata")
#load("StanCfullResults2.Rdata")
#load("StanCfullResults3.Rdata")
#load("ScoresFull.Rdata")
load("FullResultsTCGA.Rdata")


### INITIALIZE

# # create list of unique gene names and gene regions at each site
# geneNames <- list(866091)
# geneRegions <- list(866091)
# for (i in 1:866091) {
#   if (!(i %% 86609)) {
#     print(paste(i/8660.9,"% complete"))
#   }
#   geneNames[[i]] <- unique(unlist(strsplit(FullAnnotation$UCSC_RefGene_Name[i],split=";")))
#   worklist <- unique(unlist(strsplit(FullAnnotation$UCSC_RefGene_Group[i],split=";")))
#   worklist[worklist=="3'UTR"] <- "3UTR"
#   worklist[worklist=="5'UTR"] <- "5UTR"
#   worklist[worklist=="5URT"] <- "5UTR"
#   if (length(worklist)==0) { worklist <- "blank"}
#   geneRegions[[i]] <- worklist
# }
# save(geneNames, geneRegions, file = "GeneInfo.Rdata")
load("GeneInfo.Rdata")

# make vector of unique genes and regions
uniqueGenes <- unique(unlist(siteGenes))
regionTypes <- unique(unlist(siteRegions))
#uniqueGenes <- unique(unlist(geneNames))
#regionTypes <- unique(unlist(geneRegions))

# make dataframes for median and mean values of each parameter
dfMedians <- cbind.data.frame(mu_Cfull$p50, betaT_Cfull$p50, sigmaE_Cfull$p50, sigmaP_Cfull$p50, sigmaPT_Cfull$p50, sigmaT_Cfull$p50)
colnames(dfMedians) <- c("mu","betaT","sigmaE","sigmaP","sigmaPT","sigmaT")
dfMeans <- cbind.data.frame(mu_Cfull$mean, betaT_Cfull$mean, sigmaE_Cfull$mean, sigmaP_Cfull$mean, sigmaPT_Cfull$mean, sigmaT_Cfull$mean)
colnames(dfMeans) <- c("mu","betaT","sigmaE","sigmaP","sigmaPT","sigmaT")

# calculate log ratios: P/T, P/PT, PT/T
logRatios <- cbind.data.frame(log(dfMedians$sigmaP/dfMedians$sigmaT), log(dfMedians$sigmaP/dfMedians$sigmaPT), log(dfMedians$sigmaPT/dfMedians$sigmaT))
colnames(logRatios) <- c("logP/T","logP/PT","logPT/T")
logRatioMeans <- colMeans(logRatios)
# insert means by region

# create histogram plots for each log ratio, to reference later
pPT <- ggplot(logRatios, aes(x=`logP/T`)) + geom_histogram(bins = 100) + ggtitle("Histogram of logP/T") + xlim(c(-3,3)) + ylim(c(0,31000)) + geom_vline(xintercept=logRatioMeans[1], linetype="dotted")
pPPT <- ggplot(logRatios, aes(x=`logP/PT`)) + geom_histogram(bins = 100) + ggtitle("Histogram of logP/PT") + xlim(c(-3,3)) + ylim(c(0,31000)) + geom_vline(xintercept=logRatioMeans[2], linetype="dotted")
pPTT <- ggplot(logRatios, aes(x=`logPT/T`)) + geom_histogram(bins = 100) + ggtitle("Histogram of logPT/T") + xlim(c(-3,3)) + ylim(c(0,31000)) + geom_vline(xintercept=logRatioMeans[3], linetype="dotted")


### FUNCTIONS

inf2NA <- function(x) { x[is.infinite(x)] <- NA; x }

# find indices corresponding to geneString in the list geneNames
geneInd <- function (geneStr) {
  work <- grep(geneStr, geneNames)
  work2 <- sapply(work, function(i) { geneStr %in% geneNames[[i]] } )
  return(work[work2])
}
# geneInd <- function (geneStr) {
#   return(grep(paste("\\b",geneStr,"\\b",sep=""), geneNames))
# }

# plot each site along all three log ratio plots
plotLogRatios <- function (geneStr) {
  gInd <- geneInd(geneStr)
  for (i in colnames(logRatios)){
    if (i=="logP/T") { p1<-pPT }
    if (i=="logP/PT") { p1<-pPPT }
    if (i=="logPT/T") { p1<-pPTT }
    p1 <- p1 + annotate("text", x=logRatios[gInd,i], y = rep(20000,length(gInd)), label="x") + annotate("text", x=0, y=25000, label=(paste(geneStr,"mean:",format(mean(logRatios[gInd,i]),digits=5))))
    print(p1)
  }
}

# returns the functional regions of each index given
plotLogRatiosByRegions <- function (geneStr) {
  gInd <- geneInd(geneStr)
  for (i in "logP/T"){
    if (i=="logP/T") { p1<-pPT }
    if (i=="logP/PT") { p1<-pPPT }
    if (i=="logPT/T") { p1<-pPTT }
    for (j in 1:length(regionTypes)) { # go through each of the region types
      rInd <- integer(0)
      for (k in 1:length(gInd)) { # check for region j at each gInd index k
        if (regionTypes[j] %in% geneRegions[[gInd[k]]]) {
          rInd <- append(rInd,gInd[k])
        }
      }
      col1 <- c("red","orange","gold","green","blue","tomato","sienna","darkviolet")[j]
      p1 <- p1 + annotate("text", x=logRatios[rInd,i], y = rep(j*3000,length(rInd)), label="x", colour=col1) + annotate("text", x=-2, y=j*3000+1000, label=(paste(regionTypes[j],"mean:",format(mean(logRatios[rInd,i]),digits=5))),colour=col1) + annotate("text", x=2.6, y=j*3000+1000, label=(length(rInd)),colour=col1)
    }
    p1 <- p1 + annotate("text", x=logRatios[gInd,i], y = rep(27000,length(gInd)), label=".") + annotate("text", x=mean(logRatios[gInd,i]), y=28000, label=(paste(geneStr,"mean:",format(mean(logRatios[gInd,i]),digits=5)))) + annotate("text", x=logRatioMeans[i], y=30000, label=(paste("Overall mean:",format(logRatioMeans[i],digits=5)))) + annotate("text", x=2, y=29000, label=(paste("Total sites:",length(gInd))))

    print(p1)
  }
}
plotLogRatiosByRegions("GAPDH")

# score gene
scoreGene <- function(geneStr) {
  gInd <- geneInd(geneStr)
  score <- logRatioMeans*0
  for (i in colnames(logRatios)){
    score[i] <- sum(logRatios[gInd,i]-logRatioMeans[i]) / sd(logRatios[gInd,i])
  }
  return(score)
}

# score gene by functional region, including overall score
scoreGeneByRegion <- function(geneStr) {
  gInd <- geneInd(geneStr)
  score <- numeric(27)
  for (i in 1:3){
    score[i] <- (mean(logRatios[gInd,i])-logRatioMeans[i]) / sd(logRatios[gInd,i])
  }
  for (j in 1:length(regionTypes)) { # go through each of the region types
    # collect indices of region j
    rInd <- integer(0)
    for (k in 1:length(gInd)) { # check for region j at each gInd index k
      if (regionTypes[j] %in% geneRegions[[gInd[k]]]) {
        rInd <- append(rInd,gInd[k])
      }
    }
    # go through each of 3 logratios for region j
    for (i in 1:3){
      score[i+j*3] <- (mean(logRatios[rInd,i])-logRatioMeans[i]) / sd(logRatios[rInd,i])
    }
  }
  return(score)
}


### CODE

## Plots

# plot sigmas: E, P, PT, T
df1 <- dfMedians[c("sigmaE","sigmaP","sigmaPT","sigmaT")]
gg <- melt(df1)
ggplot(gg, aes(x=value, fill=variable)) +
  geom_histogram(binwidth=0.01)+xlim(c(0,2))+
  facet_grid(variable~.)

# plot log ratios
gg <- melt(logRatios)
ggplot(gg, aes(x=value, fill=variable)) +
  geom_histogram(binwidth=0.05)+
  facet_grid(variable~.)


## Check individual genes with known conservation

plotLogRatiosByRegions("APC")
plotLogRatiosByRegions("TP53")
plotLogRatiosByRegions("TTN")
plotLogRatiosByRegions("B2M")
plotLogRatiosByRegions("HLA-A")
plotLogRatiosByRegions("HLA-B")

# new genes from literature
plotLogRatiosByRegions("KRAS") # secondary mutation in CRC ***slight shift right
plotLogRatiosByRegions("PIK3CA") # oncogene ***slight shift right
plotLogRatiosByRegions("SMAD4") # colorectal tumor
plotLogRatiosByRegions("ARID1A") # TSG, colorectal, clear-cell ovarian cancer
plotLogRatiosByRegions("ATM") # TSG
plotLogRatiosByRegions("SOX9") # also related
plotLogRatiosByRegions("FAM123B") # also related ***large shift right
plotLogRatiosByRegions("MLH1") # hypermutated, mismatch repair
plotLogRatiosByRegions("MLH3") # hypermutated, mismatch repair
plotLogRatiosByRegions("MSH2") # hypermutated, mismatch repair
plotLogRatiosByRegions("MSH3") # hypermutated, mismatch repair
plotLogRatiosByRegions("MSH6") # hypermutated, mismatch repair
plotLogRatiosByRegions("PMS2") # hypermutated, mismatch repair ***slight shift right
plotLogRatiosByRegions("POLE") # hypermutated, polymerase e
plotLogRatiosByRegions("ERBB2") # recurrent copy number amplification
plotLogRatiosByRegions("IGF2") # amplified ***slight shift left
plotLogRatiosByRegions("IDH1") # oncogene ***slight shift right
plotLogRatiosByRegions("IDH2") # oncogene
plotLogRatiosByRegions("RB1") # TSG
plotLogRatiosByRegions("VHL") # TSG ***slight shift right
plotLogRatiosByRegions("NOTCH1") # oncogene/TSG ***slight shift right?
plotLogRatiosByRegions("MLL2") # medulloblastoma, therefore not related to CRC
plotLogRatiosByRegions("MLL3") # medulloblastoma, therefore not related to CRC
plotLogRatiosByRegions("SF3B1") # mRNA splicing factors, general stress
plotLogRatiosByRegions("U2AF1") # mRNA splicing factors, general stress ***slight right
plotLogRatiosByRegions("ATRX") # ALT (telomere lengthening) factors ***definite shift right
plotLogRatiosByRegions("DAXX") # ALT factors *** minimal shift right
# COSMIC genes (46)
plotLogRatiosByRegions("MDM2") # ubiquitin ligase that degrades p53, slight shift right

#random genes
plotLogRatiosByRegions("ADAM21") #
plotLogRatiosByRegions("WBSCR16") #
plotLogRatiosByRegions("TMEM213") #
plotLogRatiosByRegions("RGS18") #
plotLogRatiosByRegions("TAOK3") #
plotLogRatiosByRegions("TSPAN17") #
plotLogRatiosByRegions("COX16") #
plotLogRatiosByRegions("PTGS2") # cox2

plotLogRatiosByRegions("EN1") #

geneCount <- data.frame(table(unlist(geneNames)))
singleGenes <- geneCount$Var1[which(geneCount$Freq == 1)]

# take mean of all 3 scores, then order and remove NAs for a ranked list
geneScoresMean <- data.frame(meanscore=rowMeans(geneScores3))
rownames(geneScoresMean) <- rownames(geneScores3)
genesRanked <- geneScoresMean[order(geneScoresMean,decreasing=TRUE,na.last=FALSE),,drop=FALSE]
genesRanked2 <- genesRanked[-(1:907),,drop=FALSE]
genesRanked2$Freq <- geneCount$Freq[match(rownames(genesRanked2),geneCount$Var1)]
match(c("APC","TP53","TTN","B2M","HLA-A","HLA-B"), rownames(genesRanked2))

COSMICgenes <- c("AKT1","APC","AXIN1","AXIN2","B2M","BCL9L","BRAF","C2orf44","CTNNB1","CUX1","EIF3E","EP300","FBXW7","GRIN2A","HIF1A","KRAS","MAP2K1","MAP2K4","MDM2","MLH1","MSH2","MSH6","MUTYH","PIK3CA","PIK3R1","PMS1","PMS2","POLE","PTPRK","PTPRT","QKI","RAD21","RSPO2","RSPO3","SALL4","SFRP4","SMAD2","SMAD3","SMAD4","SRC","TBL1XR1","TCF7L2","TGFBR2","TP53","UBR5","VTI1A")
rank <- match(COSMICgenes, rownames(genesRanked2))
cbind(COSMICgenes,rank)

# analyze scores with PCA
# first 3 space (overall genes)
print(paste("Sites NA:",sum(is.na(geneScores3[,1]))))
geneScoresClear <- geneScores3[!is.na(geneScores3[,1]),]
PCAscores3 <- prcomp(geneScoresClear)
plot(PCAscores3$x[,1],PCAscores3$x[,2]) # make a scatterplot
plot(PCAscores3$x[,1],PCAscores3$x[,2], xlim=c(-4,4), ylim=c(-4,4)) # make a scatterplot
chooseRows <- (rownames(geneScoresClear) %in% c("APC","TP53","TTN","B2M","HLA-A","HLA-B"))
text(PCAscores3$x[chooseRows,1],PCAscores3$x[chooseRows,2], rownames(geneScoresClear)[chooseRows], cex=0.7, pos=4, col="red") # add labels

# 3 space again with eliminating genes with <5 sites
geneScoresClear <- geneScores3[!is.na(geneScores3[,1]),]
geneScoresClear$Freq <- geneCount$Freq[match(rownames(geneScoresClear),geneCount$Var1)]
geneScoresClear <- geneScoresClear[(geneScoresClear$Freq>5),1:3]
PCAscores3 <- prcomp(geneScoresClear)
plot(PCAscores3$x[,1],PCAscores3$x[,2]) # make a scatterplot
plot(PCAscores3$x[,1],PCAscores3$x[,2], xlim=c(-4,4), ylim=c(-4,4)) # make a scatterplot
chooseRows <- (rownames(geneScoresClear) %in% c("APC","TP53","TTN","B2M","HLA-A","HLA-B"))
text(PCAscores3$x[chooseRows,1],PCAscores3$x[chooseRows,2], rownames(geneScoresClear)[chooseRows], cex=0.7, pos=4, col="red") # add labels


# now 24 space (by functional region)      XXXX missing data probz
print(paste("Sites NA:",sum(is.na(geneScores24[,]))))
geneScoresClear <- geneScoresFull[!is.na(geneScores24[,1]),]
PCAscores24 <- prcomp(geneScores24)
geneScores24 <- replace(geneScores24, is.na(geneScores24), NA)
C <- cov(geneScores24,use="all.obs")
plot(PCAscores24$x[,1],PCAscores24$x[,2]) # make a scatterplot


hist(geneScores3$`P/T all`,3000, xlim = c(-5,5))
sum(geneScores3$`P/T all` > 0, na.rm = TRUE)/length(geneScores3$`P/T all`)
hist(geneScores3$`P/PT all`,5000, xlim = c(-5,5))
sum(geneScores3$`P/PT all` > 0, na.rm = TRUE)/length(geneScores3$`P/PT all`)
hist(geneScores3$`PT/T all`,8000, xlim = c(-5,5))
sum(geneScores3$`PT/T all` > 0, na.rm = TRUE)/length(geneScores3$`PT/T all`)



### Correlation Plots (Darryl's Suggestion)
library(gdata)
# Read in table 3 (containing genes and expressions)
tcga <- read.xls("2011-11-14592C-Sup Table 3.xls")

work <- data.frame(gene=tcga$gene,exp=tcga$median.RPKM,score=rep(NA,dim(tcga)[1]))

for (i in 1:dim(work)[1]) {
  indice <- geneInd(work$gene[i])
  work$score[i] <- mean(logRatios$`logP/T`[indice])
}

plot <- apply(is.na(work),1,sum)
ggplot(work,aes(x=log(work[,2]),y=work[,3])) + geom_point(alpha = 0.5) + labs(x="Average Expression", y="Mean Log Ratio", title="Correlation Plot")

temp <- data.frame(log(work$exp),work$score)

temp <- temp[(temp[,1] != -Inf) & (!is.na(temp[,1])) & (!is.na(temp[,2])),]
ggplot(temp,aes(x=temp[,1],y=temp[,2])) + geom_point(alpha = 0.5) + labs(x="Average Expression (log scale)", y="Mean Log Ratio", title="Correlation of Average Expression and Score")

cor(temp[,1],temp[,2])
