# script for extracting variance from Stan complex model (model3)
# input taskID corresponding to 10000 site chunks starting at site (N-1)*10000+1
# Now incorporating parallel processing by RStan methods

library(doParallel)
library(rstan)
library(gtools)

rstan_options(auto_write = FALSE)
print(paste("Cores: ", detectCores()))

# load data
load("myFA.Rdata")

args <- commandArgs(trailingOnly = TRUE)
print(args)
N <- as.numeric(args[1])
print(paste("Task #: ", N))
out.file <- args[2]
log.file <- paste("logTask",N,".txt",sep="")
writeLines(c(""), log.file)

# Choose 10000 sites by N, and select chunk within FullAnnotation, then remove FA
siteInds <- (1:10000) + (N - 1) * 10000
if (N > 86) {
  siteInds <- siteInds[1:6836]
}
nsites <- length(siteInds)

print(paste("Making chunk of nsites =", nsites))
chunk <- FullAnnotation[siteInds,]
rm(FullAnnotation)
gc()

# must add to prevent crash
checkFile <- "model3.rds"


### FXNS ###
# Function used to read in data from each site
# changed to now use a "chunk" instead of FullAnnotation, for speed and memory
# requires index 1-10000 instead of siteIndex
site <- function (ind) {
  # extract CpG site xx to start
  temp <- chunk[ind,]
  
  # here we use all patient samples, excluding glands
  indices <- c(9:13, 15:39, 46, 47, 57, 58, 72:75)
  patientLabel <- substr(colnames(temp[indices]), 1, 1)
  patientLabel[10:12] <- "K*"
  sideLabel <- substr(colnames(temp[indices]), 2, 2)
  tissueLabel <- sideLabel
  tissueLabel[tissueLabel %in% c("A", "B")] <-
    "T"   #replace A and B with T
  tumorIndicator <- 1 * (tissueLabel == "T")
  
  work <-
    data.frame(
      beta = logit(t(temp[indices])),
      patient = patientLabel,
      tissue = tissueLabel,
      side = sideLabel,
      tInd = tumorIndicator
    )
  
  return(work)
}

# Function to build stan model for each site
stanfit3 <- function (dataset) {
  stanDat <- list(
    pID = as.integer(factor(dataset$patient)),
    tInd = dataset$tInd,
    N = nrow(dataset),
    P = nlevels(dataset$patient),
    y = dataset[, 1]
  )
  
  # Using Model 3: add intra-tumoral variances
  # have to check for model3.rds - crashes Stan
  if (file.exists(checkFile)) {
    file.remove(checkFile)
    print("caught one!") 
  }
  stanFit3 <-
    stan(
      #file = "model3.stan",
      fit = emptyFit,
      data = stanDat,
      control = list(adapt_delta = 0.999),
      refresh = 0
    )
  
  return(stanFit3 = stanFit3)
}

# Function for combining parallel runs of each fixef/sigma
comb <- function(x, ...) {
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}


### CODE ###

print("Generating data.frames")

betaT_C <-
  data.frame(
    mean = numeric(nsites),
    p2.5 = numeric(nsites),
    p25 = numeric(nsites),
    p50 = numeric(nsites),
    p75 = numeric(nsites),
    p97.5 = numeric(nsites)
  )

mu_C <-
  data.frame(
    mean = numeric(nsites),
    p2.5 = numeric(nsites),
    p25 = numeric(nsites),
    p50 = numeric(nsites),
    p75 = numeric(nsites),
    p97.5 = numeric(nsites)
  )

sigmaE_C <-
  data.frame(
    mean = numeric(nsites),
    p2.5 = numeric(nsites),
    p25 = numeric(nsites),
    p50 = numeric(nsites),
    p75 = numeric(nsites),
    p97.5 = numeric(nsites)
  )

sigmaP_C <-
  data.frame(
    mean = numeric(nsites),
    p2.5 = numeric(nsites),
    p25 = numeric(nsites),
    p50 = numeric(nsites),
    p75 = numeric(nsites),
    p97.5 = numeric(nsites)
  )

sigmaT_C <-
  data.frame(
    mean = numeric(nsites),
    p2.5 = numeric(nsites),
    p25 = numeric(nsites),
    p50 = numeric(nsites),
    p75 = numeric(nsites),
    p97.5 = numeric(nsites)
  )

sigmaPT_C <-
  data.frame(
    mean = numeric(nsites),
    p2.5 = numeric(nsites),
    p25 = numeric(nsites),
    p50 = numeric(nsites),
    p75 = numeric(nsites),
    p97.5 = numeric(nsites)
  )

# inds for extracting data
mInd <- c(1, 4:8)

# create empty fit to use same model compilation throughout
data1<-site(1)
stanDat1 <- list(
  pID = as.integer(factor(data1$patient)),
  tInd = data1$tInd,
  N = nrow(data1),
  P = nlevels(data1$patient),
  y = data1[, 1]
)
emptyFit <- stan(file="model3.stan", data = stanDat1, chains = 0)

# run in parallel via doParallel
cl <- makeCluster(detectCores())
registerDoParallel(cl)
print(paste("Cores registered:",getDoParWorkers()))
print(paste("Backend type:",getDoParName()))
print("Starting foreach loop")
ptm <- proc.time()
parData <- foreach(i = iter(1:nsites), .combine = 'comb', .multicombine = TRUE, .packages=c("gtools","rstan")) %dopar% {
  sink(log.file, append=TRUE)
  print(paste("site:", i,"/",nsites))
  sink()

  data <- site(i)
  stanFit <- stanfit3(data)
  fitSumm <- summary(stanFit)$summary[71:76, mInd]
  rm(data,stanFit)
  gc()
  
  split(fitSumm, row(fitSumm))
}
proc.time() - ptm
stopCluster(cl)

betaT_C[,] <- parData$'1'
mu_C[,] <- parData$'2'
sigmaE_C[,] <- parData$'3'
sigmaP_C[,] <- parData$'4'
sigmaPT_C[,] <- parData$'5'
sigmaT_C[,] <- parData$'6'
rm(parData)
gc()

# # parallel via Rstan
# options(mc.cores = detectCores())
# ptm <- proc.time()
# for (i in (1:nsites)) {
#   sink(log.file, append=TRUE)
#   print(paste("site:", i))
#   sink()
#   
#   data <- site(siteInds[i])
#   stanFit <- stanfit3(data)
#   fitSumm <- summary(stanFit)$summary[71:76,mInd]
#   rm(data,stanFit)
#   gc()
#   
#   betaT_C[i,] <- fitSumm[1,]
#   mu_C[i,] <- fitSumm[2,]
#   sigmaE_C[i,] <- fitSumm[3,]
#   sigmaP_C[i,] <- fitSumm[4,]
#   sigmaPT_C[i,] <- fitSumm[5,]
#   sigmaT_C[i,] <- fitSumm[6,]
# }
# proc.time() - ptm

print(paste("Completed run, now saving"))
# save data
save(mu_C, betaT_C, sigmaP_C, sigmaT_C,
     sigmaPT_C, sigmaE_C, file = out.file)
