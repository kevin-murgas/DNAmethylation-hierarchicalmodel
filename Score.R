########################## Scoring ################################

library(doParallel)

### Initialize

print(paste("Cores: ", detectCores()))

# here set the working directory that points to the data folder
# e.g. the folder with annotated data saved as "myFA.Rdata"
# please just comment out the directories already here, so each
# user can uncomment the setwd for them as we move code around

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# Yanlin's working directory
setwd("D:/DataPlus2017/Data")


### LOAD DATA

# load fully annotated data (saved from LoadDataAndQC.R)
load("myFA.Rdata")
load("GeneInfo.Rdata")
load("StanCfullResults3.Rdata")

N<-1

# make vector of unique genes and regions
uniqueGenes <- unique(unlist(geneNames))
regionTypes <- unique(unlist(geneRegions))

# make dataframes for median values of each parameter
dfMedians <- cbind.data.frame(mu_Cfull3$p50, betaT_Cfull3$p50, sigmaE_Cfull3$p50, sigmaP_Cfull3$p50, sigmaPT_Cfull3$p50, sigmaT_Cfull3$p50)
colnames(dfMedians) <- c("mu","betaT","sigmaE","sigmaP","sigmaPT","sigmaT")

# calculate log ratios: P/T, P/PT, PT/T
logRatios <- cbind.data.frame(log(dfMedians$sigmaP/dfMedians$sigmaT), log(dfMedians$sigmaP/dfMedians$sigmaPT), log(dfMedians$sigmaPT/dfMedians$sigmaT))
colnames(logRatios) <- c("logP/T","logP/PT","logPT/T")
logRatioMeans <- colMeans(logRatios)
# insert means by region

# Choose 1000 genes from uniqueGenes by N, and select geneChunk within gene list
geneIDs <- (1:1000) + (N - 1) * 1000
if (N > 27) {
  geneIDs <- geneIDs[1:383]
}
geneIDs <- 1:20
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

# score gene
# five methods: 
# 1. mean distance from overall mean, averaged over all sites
# 2. mean distance from overall mean divided by standard deviation of gene
# 3. average PTprob of all sites in gene
# 4. fraction of sites in gene with PTprob > 0.95 (binary)
# 5. Higher Criticism Statistic
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

log.file <- paste0("logTask",N,".txt")
writeLines(c(""), log.file)

## score all genes
# run in parallel via doParallel
cl <- makeCluster(detectCores()-2)
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
rownames(geneScores) <- geneChunk
colnames(geneScores) <- c(1:5, "siteCount")

cbind(geneScores, geneCount$Freq[1:20])
plot(geneScores[,4], geneScores[,5])

print(paste("Completed run, now saving"))
# save data
save(geneScores, file = out.file)
