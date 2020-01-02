# load StanCRun data
# produce full size data frames for mu,betaT,sigmas (P,PT,T,E)
# read in individual save files containing 10000 site runs
# store 10000 site chunks in corresponding positions in full data frame
# repeat for each completed run

library(ggplot2)
library(reshape2)

### LOADING FOR STAN MODEL RESULTS

# Set working directory to normal data wd then folder ResultsC
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data/NewResultsC")

betaT_Cfull3 <-
  data.frame(
    mean = numeric(866836),
    p2.5 = numeric(866836),
    p25 = numeric(866836),
    p50 = numeric(866836),
    p75 = numeric(866836),
    p97.5 = numeric(866836)
  )

mu_Cfull3 <- betaT_Cfull3
sigmaE_Cfull3 <- betaT_Cfull3
sigmaP_Cfull3 <- betaT_Cfull3
sigmaPT_Cfull3 <- betaT_Cfull3
sigmaT_Cfull3 <- betaT_Cfull3
PTprob3 <- numeric(866836)

nsites <- 10000
filesPresent <- (1:87)
fileNames <- dir(pattern = "StanCParResults_*")
allInds <- numeric()
for (i in filesPresent) {
  print(paste("Chunk:",i))
  load(fileNames[fileNames == paste("StanCParResults_",i,".Rdata", sep="")])
  inds <- (1:nsites) + (i - 1) * nsites
  if (i==87) {
    inds <- inds[1:6836]
  }
  allInds <- append(allInds, inds)
  betaT_Cfull3[inds, ] <- betaT_C
  mu_Cfull3[inds, ] <- mu_C
  sigmaE_Cfull3[inds, ] <- sigmaE_C
  sigmaP_Cfull3[inds, ] <- sigmaP_C
  sigmaPT_Cfull3[inds, ] <- sigmaPT_C
  sigmaT_Cfull3[inds, ] <- sigmaT_C
  PTprob3[inds] <- PTprob
}
save(mu_Cfull3, betaT_Cfull3, sigmaP_Cfull3, sigmaT_Cfull3,
     sigmaPT_Cfull3, sigmaE_Cfull3, PTprob3, file = "StanCfullResults3.Rdata")


### LOADING FOR GENE SCORING RESULTS
# used to be for 3- and 24-space
# now for 5 score methods on logPTratio only

# Set working directory to normal data wd then folder GeneScores2
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data/GeneScores2")

geneScores5 <- as.data.frame(matrix(ncol = 6, nrow = 27383))
#geneScores27 <- as.data.frame(matrix(ncol = 27, nrow = 27383))
#colnames(geneScores27) <- paste(rep(c("P/T","P/PT","PT/T"),8),c(rep("all",3),rep(regionTypes[1],3),rep(regionTypes[2],3),rep(regionTypes[3],3),rep(regionTypes[4],3),rep(regionTypes[5],3),rep(regionTypes[6],3),rep(regionTypes[7],3),rep(regionTypes[8],3)))
#colnames(geneScores3) <- c("logP/T","logP/PT","logPT/T")

nGenes <- 1000
filesPresent <- (1:28)
fileNames <- dir(pattern = "ScoreResults_*")
allInds <- numeric()
for (i in filesPresent) {
  print(paste("Chunk:",i))
  load(fileNames[fileNames == paste0("ScoreResults_",i,".Rdata")])
  inds <- (1:nGenes) + (i - 1) * nGenes
  if (i==28) {
    inds <- inds[1:383]
  }
  allInds <- append(allInds, inds)
  geneScores5[inds,] <- geneScores
  rownames(geneScores5)[inds] <- rownames(geneScores)
}
rownames(geneScores5) <- uniqueGenes # overwrite due to improper geneInfo file
colnames(geneScores5) <- c("MeanDist 1", "MeanDist/SD 2", "MeanPTProb 3", "FracBinProb 4", "HigherCrit 5", "Site Count")
# five methods: 
# 1. mean distance from overall mean, averaged over all sites
# 2. mean distance from overall mean divided by standard deviation of gene
# 3. average PTprob of all sites in gene
# 4. fraction of sites in gene with PTprob > 0.95 (binary)
# 5. Higher Criticism Statistic

# since fixed but the saved results include dates instead of numbers for number-name genes

save(geneScores5, file = "ScoresFull5.Rdata")
write.csv(geneScores5, file = "GeneScores5Methods.csv")



### compare modelC runs 1 and 2
# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")
load("StanCfullResults.Rdata")
load("StanCfullResults2.Rdata")

mu1 <- mu_Cfull$p50
mu2 <- mu_Cfull2$p50
cor(mu1,mu2)
mudiffs <- mu1-mu2
muchange <- mudiffs/mu1
summary(mudiffs)
summary(muchange)
hist(mudiffs,100)
hist(muchange,100)

betaT1 <- betaT_Cfull$p50
betaT2 <- betaT_Cfull2$p50
cor(betaT1,betaT2)
betaTdiffs <- betaT1-betaT2
betaTchange <- betaTdiffs/betaT1
summary(betaTdiffs)
summary(betaTchange)
hist(betaTdiffs,100)
hist(betaTchange,100)

sigmaE1 <- sigmaE_Cfull$p50
sigmaE2 <- sigmaE_Cfull2$p50
cor(sigmaE1,sigmaE2)
sigmaEdiffs <- sigmaE1-sigmaE2
sigmaEchange <- sigmaEdiffs/sigmaE1
summary(sigmaEdiffs)
summary(sigmaEchange)
hist(sigmaEdiffs,100)
hist(sigmaEchange,100)

sigmaP1 <- sigmaP_Cfull$p50
sigmaP2 <- sigmaP_Cfull2$p50
cor(sigmaP1,sigmaP2)
sigmaPdiffs <- sigmaP1-sigmaP2
sigmaPchange <- sigmaPdiffs/sigmaP1
summary(sigmaPdiffs)
summary(sigmaPchange)
hist(sigmaPdiffs,100)
hist(sigmaPchange,100)

sigmaPT1 <- sigmaPT_Cfull$p50
sigmaPT2 <- sigmaPT_Cfull2$p50
cor(sigmaPT1,sigmaPT2)
sigmaPTdiffs <- sigmaPT1-sigmaPT2
sigmaPTchange <- sigmaPTdiffs/sigmaPT1
summary(sigmaPTdiffs)
summary(sigmaPTchange)
hist(sigmaPTdiffs,100)
hist(sigmaPTchange,100)

sigmaT1 <- sigmaT_Cfull$p50
sigmaT2 <- sigmaT_Cfull2$p50
cor(sigmaT1,sigmaT2)
sigmaTdiffs <- sigmaT1-sigmaT2
sigmaTchange <- sigmaTdiffs/sigmaT1
summary(sigmaTdiffs)
quantile(sigmaTdiffs, c(0,0.05,0.25,0.5,0.75,0.95,1))
summary(sigmaTchange)
format(quantile(sigmaTchange, c(0,0.05,0.25,0.5,0.75,0.95,1)),digits=4)
hist(sigmaTdiffs,100)
hist(sigmaTchange,100)


