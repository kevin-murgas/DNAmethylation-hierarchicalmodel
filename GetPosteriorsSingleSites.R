# script for extracting variance from Stan complex model (model3)
# input taskID corresponding to 10000 site chunks starting at site (N-1)*10000+1
# Now incorporating parallel processing by RStan methods

library(doParallel)
library(rstan)
library(gtools)
library(ggplot2)
library(bayesplot)


## INITIALIZE
rstan_options(auto_write = FALSE)
print(paste("Cores: ", detectCores()))

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# Yanlin's working directory
#setwd("D:/DataPlus2017/Data")

# load data
load("myFA.Rdata")

#args <- commandArgs(trailingOnly = TRUE)
#print(args)
#N <- as.numeric(args[1])
N <- 1
print(paste("Task #: ", N))
#out.file <- args[2]
log.file <- paste("logTask",N,".txt",sep="")
writeLines(c(""), log.file)

# randomly select 5 sites in 1-866836
siteInds <- sample(1:866836,5)
nsites <- length(siteInds)

# select chunk within FullAnnotation, then remove FA
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
      warmup = 1000,
      iter = 4000,
      chains = 4,
      thin = 1,
      data = stanDat,
      control = list(adapt_delta = 0.999),
      refresh = 0
    )
  
  return(stanFit3 = stanFit3)
}


### CODE ###

# check for model3.rds
if (file.exists(checkFile)) {
  file.remove(checkFile)
  print("caught one!") 
}

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

# # run in parallel via doParallel
# # want to collect data of random effect coefficients (bcd), fixed effects (betaT,mu), and sigmas
# cl <- makeCluster(detectCores()-3)
# registerDoParallel(cl)
# print(paste("Cores registered:",getDoParWorkers()))
# print(paste("Backend type:",getDoParName()))
# print("Starting foreach loop")
# ptm <- proc.time()
# parData <- foreach(i = iter(1:nsites), .combine = 'comb', .multicombine = TRUE, .packages=c("gtools","rstan")) %dopar% {
#   sink(log.file, append=TRUE)
#   print(paste("site:", i,"/",nsites))
#   sink()
# 
#   data <- site(i)
#   stanFit <- stanfit3(data)
#   fitSumm <- summary(stanFit)$summary[71:76, mInd] # pulls betaT, mu, sigmaE,P,PT,T
# 
#   posterior <- as.matrix(stanFit,pars=c("sigma_p","sigma_t"))
#   PTratio <- log(posterior[,1]/posterior[,2])
#   prob <- sum(PTratio > 0)/length(PTratio)
#   fitSumm <- rbind(fitSumm,rep(prob,length(mInd)))
# 
#   rm(data,stanFit)
#   gc()
#   
#   split(fitSumm, row(fitSumm))
# }
# proc.time() - ptm
# stopCluster(cl)
# 
# betaT_C[,] <- parData$'1'
# mu_C[,] <- parData$'2'
# sigmaE_C[,] <- parData$'3'
# sigmaP_C[,] <- parData$'4'
# sigmaPT_C[,] <- parData$'5'
# sigmaT_C[,] <- parData$'6'
# PTprob <- parData$'7'[,1]
# #rm(parData)
# gc()

# run in serial
# want to collect data of random effect coefficients (bcd), fixed effects (betaT,mu), and sigmas
# going to make it run the site and present the data within the same for loop, for data efficiency
ptm <- proc.time()
for(i in 1:nsites) {
  sink(log.file, append=TRUE)
  print(paste("site:", i,"/",nsites))
  sink()
  
  data <- site(i)
  colnames(data)[1] <- "beta"
  stanFit <- stanfit3(data)
  
  # plot beta and fits at this site for all patients, using mean estimates for b and c coefficients
  p1 <- ggplot() + geom_point(data = data, aes(x=tInd, y=beta, colour=patient)) + ggtitle(paste("Site",siteInds[i],": Beta Fits by Patient"))
  pats <- as.integer(factor(data$patient))
  npats <- max(pats)
  patcoefs <- data.frame(pat = rep(1:npats), b = summary(stanFit)$summary[1:npats], c = summary(stanFit)$summary[(1:npats)+npats])
  betaT <- summary(stanFit)$summary[71,1]
  mu <- summary(stanFit)$summary[72,1]
  linedf <- data.frame(pat = rep(levels(data$patient)[1:npats],2), tInd = c(rep(0,npats),rep(1,npats)), est = c(patcoefs$b+mu, (patcoefs$b+patcoefs$c+mu+betaT)))
  p1 <- p1 + geom_line(data=linedf, aes(x=tInd, y=est, colour=pat))
  
  # plot posteriors for this site
  # variables: mu, betaT, sigmaP, sigmaT, sigmaPT, sigmaE
  posterior <- as.array(stanFit)[,,71:76]
  p2 <- mcmc_areas(
    posterior, 
    pars = c("mu", "betaT", "sigma_e", "sigma_p", "sigma_t","sigma_pt"),
    prob = 0.8, # 80% intervals
    prob_outer = 0.95, 
    point_est = "mean"
  ) + ggtitle(paste("Site",siteInds[i],": Fixed Effect and Sigma Posteriors"))
  
  # also plot PTratio?
  posteriorPT <- as.array(stanFit,pars=c("sigma_p","sigma_t"))
  PTratio <- log(posteriorPT[,,1, drop=FALSE]/posteriorPT[,,2, drop=FALSE])
  prob <- sum(PTratio > 0)/length(PTratio)
  p3 <- mcmc_areas(
    PTratio,
    prob = 0.8, # 80% intervals
    prob_outer = 0.95,
    point_est = "mean"
  ) + ggtitle(paste("Site",siteInds[i],": LogPTratio Posterior (PTprob = ",prob,")"))
  
  print(p1)
  print(p2)
  print(p3)
  
  rm(data,stanFit,posterior)
  gc()
}
proc.time() - ptm

print(paste("Completed run on sites:", paste(siteInds, collapse=" ")))

