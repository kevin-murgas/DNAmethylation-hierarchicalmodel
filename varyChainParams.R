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
siteInds <- sample(1:866836,1)
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
# chain 1 - default; 2000 warmup, 4000 total iter, 4 chains, thin 1
stanfit3chain1 <- function (dataset) {
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
      warmup = 2000,
      iter = 4000,
      chains = 4,
      thin = 1,
      data = stanDat,
      control = list(adapt_delta = 0.999),
      refresh = 0
    )
  
  return(stanFit3 = stanFit3)
}

# chain 2 - long chain; 2000 warmup, 10000 total iter, 2 chains, thin 1
stanfit3chain2 <- function (dataset) {
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
      warmup = 2000,
      iter = 10000,
      chains = 2,
      thin = 1,
      data = stanDat,
      control = list(adapt_delta = 0.999),
      refresh = 0
    )
  
  return(stanFit3 = stanFit3)
}

# chain 3 - short warmup; 1000 warmup, 4000 total iter, 4 chains, thin 1
stanfit3chain3 <- function (dataset) {
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

# chain 4 - many chains; 2000 warmup, 4000 total iter, 10 chains, thin 1
stanfit3chain4 <- function (dataset) {
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
      warmup = 2000,
      iter = 4000,
      chains = 10,
      thin = 1,
      data = stanDat,
      control = list(adapt_delta = 0.999),
      refresh = 0
    )
  
  return(stanFit3 = stanFit3)
}

### CODE ###

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
# want plots for each chain (title with chain number):
# sigmas E,P,PT,T posterior distributions
# sigmas E,P,PT,T trace plots- examine with warmup once to see if worth (prob not)
# log-posterior(P/T) posterior distribution and traceplot
# lp? measure of error?
# 
ptmAll <- proc.time()
for(i in 1:nsites) {
  sink(log.file, append=TRUE)
  print(paste("site:", i,"/",nsites))
  sink()
  
  data <- site(i)
  colnames(data)[1] <- "beta"
  print("begin chain 1")
  ptm <- proc.time()
  stanFitChain1 <- stanfit3chain1(data)
  print(proc.time() - ptm)
  print("end chain 1")
  print("begin chain 2")
  ptm <- proc.time()
  stanFitChain2 <- stanfit3chain2(data)
  print(proc.time() - ptm)
  print("end chain 2")
  print("begin chain 3")
  ptm <- proc.time()
  stanFitChain3 <- stanfit3chain3(data)
  print(proc.time() - ptm)
  print("end chain 3")
  print("begin chain 4")
  ptm <- proc.time()
  stanFitChain4 <- stanfit3chain4(data)
  print(proc.time() - ptm)
  print("end chain 4")

  # plot beta and fits at this site for all patients, using mean estimates for b and c coefficients
  # start with p1 just data points
  p1 <- ggplot() + geom_point(data = data, aes(x=tInd, y=beta, colour=patient)) + ggtitle(paste("Site",siteInds[i],": Beta Fits by Patient"))
  pats <- as.integer(factor(data$patient))
  npats <- max(pats)
  #ch1
  patcoefs <- data.frame(pat = rep(1:npats), b = summary(stanFitChain1)$summary[1:npats], c = summary(stanFitChain1)$summary[(1:npats)+npats])
  betaT <- summary(stanFitChain1)$summary[71,1]
  mu <- summary(stanFitChain1)$summary[72,1]
  linedf <- data.frame(pat = rep(levels(data$patient)[1:npats],2), tInd = c(rep(0,npats),rep(1,npats)), est = c(patcoefs$b+mu, (patcoefs$b+patcoefs$c+mu+betaT)))
  p1 <- p1 + geom_line(data=linedf, aes(x=tInd, y=est, colour=pat))
  #ch2
  patcoefs <- data.frame(pat = rep(1:npats), b = summary(stanFitChain2)$summary[1:npats], c = summary(stanFitChain2)$summary[(1:npats)+npats])
  betaT <- summary(stanFitChain2)$summary[71,1]
  mu <- summary(stanFitChain2)$summary[72,1]
  linedf <- data.frame(pat = rep(levels(data$patient)[1:npats],2), tInd = c(rep(0,npats),rep(1,npats)), est = c(patcoefs$b+mu, (patcoefs$b+patcoefs$c+mu+betaT)))
  p1 <- p1 + geom_line(data=linedf, aes(x=tInd, y=est, colour=pat))
  #ch3
  patcoefs <- data.frame(pat = rep(1:npats), b = summary(stanFitChain3)$summary[1:npats], c = summary(stanFitChain3)$summary[(1:npats)+npats])
  betaT <- summary(stanFitChain3)$summary[71,1]
  mu <- summary(stanFitChain3)$summary[72,1]
  linedf <- data.frame(pat = rep(levels(data$patient)[1:npats],2), tInd = c(rep(0,npats),rep(1,npats)), est = c(patcoefs$b+mu, (patcoefs$b+patcoefs$c+mu+betaT)))
  p1 <- p1 + geom_line(data=linedf, aes(x=tInd, y=est, colour=pat))
  #ch4
  patcoefs <- data.frame(pat = rep(1:npats), b = summary(stanFitChain4)$summary[1:npats], c = summary(stanFitChain4)$summary[(1:npats)+npats])
  betaT <- summary(stanFitChain4)$summary[71,1]
  mu <- summary(stanFitChain4)$summary[72,1]
  linedf <- data.frame(pat = rep(levels(data$patient)[1:npats],2), tInd = c(rep(0,npats),rep(1,npats)), est = c(patcoefs$b+mu, (patcoefs$b+patcoefs$c+mu+betaT)))
  p1 <- p1 + geom_line(data=linedf, aes(x=tInd, y=est, colour=pat))
  
  print(p1) # all chains monster plot, not actually that monstrous
  # print(p1chain1)
  # print(p1chain2)
  # print(p1chain3)
  # print(p1chain4)
  
  
  # plot trace plots and posteriors for this site for sigma_e, sigma_p, sigma_t, sigma_pt
  p1 <- traceplot(stanFitChain1, pars = c("sigma_e","sigma_p","sigma_pt","sigma_t"))
  print(p1)
  # p1 <- traceplot(stanFitChain1, pars = c("lp__")) #, inc_warmup = TRUE 
  # print(p1)
  posteriorChain1 <- as.array(stanFitChain1)[,,71:76]
  p2chain1 <- mcmc_areas(
    posteriorChain1, pars = c("sigma_e", "sigma_p", "sigma_t","sigma_pt"),
    prob = 0.8,prob_outer = 0.95,point_est = "mean"
  ) + ggtitle(paste("Site",siteInds[i],": Fixed Effect and Sigma Posteriors, Chain1"))
  #ch2
  p1 <- traceplot(stanFitChain2, pars = c("sigma_e","sigma_p","sigma_pt","sigma_t"))
  print(p1)
  posteriorChain2 <- as.array(stanFitChain2)[,,71:76]
  p2chain2 <- mcmc_areas(
    posteriorChain2, pars = c("sigma_e", "sigma_p", "sigma_t","sigma_pt"),
    prob = 0.8,prob_outer = 0.95,point_est = "mean"
  ) + ggtitle(paste("Site",siteInds[i],": Fixed Effect and Sigma Posteriors, Chain2"))
  # ch3
  p1 <- traceplot(stanFitChain3, pars = c("sigma_e","sigma_p","sigma_pt","sigma_t"))
  print(p1)
  posteriorChain3 <- as.array(stanFitChain3)[,,71:76]
  p2chain3 <- mcmc_areas(
    posteriorChain3, pars = c("sigma_e", "sigma_p", "sigma_t","sigma_pt"),
    prob = 0.8,prob_outer = 0.95,point_est = "mean"
  ) + ggtitle(paste("Site",siteInds[i],": Fixed Effect and Sigma Posteriors, Chain3"))
  # ch4
  p1 <- traceplot(stanFitChain4, pars = c("sigma_e","sigma_p","sigma_pt","sigma_t"))
  print(p1)
  posteriorChain4 <- as.array(stanFitChain4)[,,71:76]
  p2chain4 <- mcmc_areas(
    posteriorChain4, pars = c("sigma_e", "sigma_p", "sigma_t","sigma_pt"),
    prob = 0.8,prob_outer = 0.95,point_est = "mean"
  ) + ggtitle(paste("Site",siteInds[i],": Fixed Effect and Sigma Posteriors, Chain4"))
  
  print(p2chain1)
  print(p2chain2)
  print(p2chain3)
  print(p2chain4)
  
  rm(data)
  rm(stanFitChain1,posteriorChain1)
  rm(stanFitChain2,posteriorChain2)
  rm(stanFitChain3,posteriorChain3)
  rm(stanFitChain4,posteriorChain4)
  gc()
}
proc.time() - ptmAll

print(paste("Completed run on sites:", paste(siteInds, collapse=" ")))

