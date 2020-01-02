# script for calculating scores on individual genes (27383 total genes)
# using full-scale results from Stan complex model (model3), including PTprob (run 3)
# input taskID (1-28) corresponding to 1000 gene chunks starting at site (N-1)*1000+1
# incorporating parallel processing to run through all genes faster
# just need "StanCfullResults3.Rdata" and "GeneInfo" in working directory

library(doParallel)

### Initialize

print(paste("Cores: ", detectCores()))

args <- commandArgs(trailingOnly = TRUE)
print(args)
N <- as.numeric(args[1])
print(paste("Task #:", N))
out.file <- args[2]
log.file <- paste0("logTask",N,".txt")
writeLines(c(""), log.file)

# load data: Stan parameters, geneNames and geneRegions lists
load("StanCfullResults3.Rdata")
load("GeneInfo.Rdata")

# make vector of unique genes and regions
uniqueGenes <- unique(unlist(geneNames))
regionTypes <- unique(unlist(geneRegions))

# count sites at each geneName
geneCount <- data.frame(table(unlist(geneNames)))

# make dataframes for median values of each parameter
dfMedians <- cbind.data.frame(mu_Cfull$p50, betaT_Cfull$p50, sigmaE_Cfull$p50, sigmaP_Cfull$p50, sigmaPT_Cfull$p50, sigmaT_Cfull$p50)
colnames(dfMedians) <- c("mu","betaT","sigmaE","sigmaP","sigmaPT","sigmaT")

# calculate log ratios: P/T, P/PT, PT/T
logRatios <- cbind.data.frame(log(dfMedians$sigmaP/dfMedians$sigmaT), log(dfMedians$sigmaP/dfMedians$sigmaPT), log(dfMedians$sigmaPT/dfMedians$sigmaT))
colnames(logRatios) <- c("logP/T","logP/PT","logPT/T")
logRatioMeans <- colMeans(logRatios)

# Choose 1000 genes from uniqueGenes by N, and select geneChunk within gene list
geneIDs <- (1:1000) + (N - 1) * 1000
if (N > 27) {
  geneIDs <- geneIDs[1:383]
}
numGenes <- length(geneIDs)

print(paste("Making chunk of numGenes =", numGenes))
geneChunk <- uniqueGenes[geneIDs]


### Functions

# find indices corresponding to geneString in the list geneNames
geneInd <- function (geneStr) {
  work <- grep(geneStr, geneNames)
  work2 <- sapply(work, function(i) { geneStr %in% geneNames[[i]] } )
  return(work[work2])
}

# score gene by 5 methods, using logPTratio
scoreGene <- function(geneStr) {
  gInd <- geneInd(geneStr)
  score <- numeric(6)
  
  # 1. mean distance from overall mean, averaged over all sites
  score[1] <- mean(logRatios[gInd,1]-logRatioMeans[1])
  
  # 2. mean distance from overall mean divided by standard deviation of gene
  score[2] <- mean(logRatios[gInd,1]-logRatioMeans[1]) / sd(logRatios[gInd,1])
  
  # 3. average PTprob of all sites in gene
  score[3] <- mean(PTprob3[gInd])
  
  # 4. fraction of sites in gene with PTprob > 0.95 (binary)
  binprob <- (PTprob3 > 0.95)
  score[4] <- mean(binprob[gInd])
  
  # 5. Higher Criticism Statistic
  n <- length(gInd)
  sortprob <- sort(PTprob3[gInd])
  # catch any 0/1's to prevent Inf
  sortprob[sortprob==0] <- 0.001
  sortprob[sortprob==1] <- 0.999
  normz <- sqrt(n)*abs((1:n)/n - sortprob)/(sqrt(sortprob*(1-sortprob)))
  score[5] <- max(normz)
  
  score[6] <- n
  
  return(score)
}


### Code

## score all genes
# run in parallel via doParallel
cl <- makeCluster(detectCores())
registerDoParallel(cl)
print(paste("Cores registered:",getDoParWorkers()))
print(paste("Backend type:",getDoParName()))
print("Starting foreach loop")
ptm <- proc.time()
parData <- foreach(g=iter(1:numGenes), .combine=rbind) %dopar% {
  sink(log.file, append=TRUE)
  print(paste("Gene:", g,"/",numGenes))
  sink()
  scoreGene(geneChunk[g])
}
proc.time() - ptm
stopCluster(cl)

geneScores <- parData
rm(parData)
gc()
rownames(geneScores) <- geneChunk
colnames(geneScores) <- c(1:5, "siteCount")

print(paste("Completed run, now saving"))
# save data
save(geneScores, file = out.file)
