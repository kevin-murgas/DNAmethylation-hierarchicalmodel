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
load("myFA_bulkonly.Rdata")

# randomly select 5 sites in 1-866091
#siteInds <- sample(1:866091,5)
siteInds <- c(525490, 570532, 827025, 287316, 854143)
siteInds <- 124022
nsites <- length(siteInds)


### FXNS ###
# Function used to read in data from each site
site <- function (ind) {
  # extract CpG site xx to start
  temp <- FullAnnotation[ind,]
  
  # here we use all patient samples, excluding glands
  #indices <- c(10:14, 16:40,47,48,58,59,73:76,87:96) # if using full data file
  indices <- c(10:57) # if using bulk only data file
  patientLabel <- substr(colnames(temp[indices]), 1, 1)
  patientLabel[10:12] <- "K*"
  sideLabel <- substr(colnames(temp[indices]), 2, 2)
  tissueLabel <- sideLabel
  tissueLabel[!(sideLabel %in% c("N"))] <- "T"   #replace A/B/D/M with T
  tumorIndicator <- 1 * (tissueLabel == "T")
  
  work <-
    data.frame(
      M = logit(t(temp[indices])),
      patient = patientLabel,
      tissue = tissueLabel,
      side = sideLabel,
      tInd = tumorIndicator
    )
  colnames(work)[1] <- "M"
  return(work)
}

# Function to build stan model for each site
stanfit_model <- function (dataset) {
  stanDat <- list(
    N = nrow(dataset),
    P = nlevels(dataset$patient),
    y = dataset[, 1],
    tInd = dataset$tInd,
    pID = as.integer(factor(dataset$patient))
  )
  stanFit <- 
    stan(
      fit = emptyFit,
      warmup = 200,
      iter = 2000,
      chains = 4,
      data = stanDat,
      control = list(adapt_delta = 0.999),
      refresh = 0,
      seed=123
    )
  return(stanFit=stanFit)
}

# Plot raw data and fits
fitplot <- function(data,stanFit) {
  # plot beta and fits at this site for all patients, using mean estimates for b and c coefficients
  # first plot data points
  p1 <- ggplot() + geom_point(data = data, aes(x=tInd, y=beta, colour=patient))
  
  # extract coefficients from stan fit
  pats <- as.integer(factor(data$patient))
  nsamp <- length(pats)
  npats <- max(pats)
  sumTemp <- summary(stanFit)$summary
  mu <- sumTemp[1,1] # mean estimate
  betaT <- sumTemp[2,1] # mean
  sampcoefs <- data.frame(pat = pats, b = sumTemp[6+pats], c = sumTemp[6+npats+pats], d = sumTemp[6+2*npats+c(1:nsamp)])
  patcoefs <- data.frame(pat = 1:npats, b = sumTemp[6+1:npats], c = sumTemp[6+npats+1:npats])
  rm(sumTemp)
  # add individual sample estimates
  sampdf <- data.frame(pat = rep(factor(data$patient),2), patg = rep(factor(1:nsamp),2), tInd = c(rep(0,nsamp),rep(1,nsamp)), est = c(sampcoefs$b+mu, (sampcoefs$b+sampcoefs$c++sampcoefs$d+mu+betaT)))
  p1 <- p1 + geom_line(data=sampdf, aes(x=tInd, y=est, colour=pat, group=patg), linetype="dashed", alpha=0.5)
  # add patient mean estimates
  patdf <- data.frame(pat = rep(levels(data$patient),2), tInd = c(rep(0,npats),rep(1,npats)), est = c(patcoefs$b+mu, patcoefs$b+patcoefs$c+mu+betaT))
  p1 <- p1 + geom_line(data=patdf, aes(x=tInd, y=est, colour=pat))
  # add overall mean estimates as black dots
  meandf <- data.frame(tInd = c(0,1), est = c(mu,mu+betaT))
  p1 <- p1 + geom_point(data=meandf, aes(x=tInd, y=est), size=3, alpha=0.7) + geom_line(data=meandf, aes(x=tInd, y=est), size=2, alpha=0.5)
  
  return(p1)
}

# Plot posterior distribution
postplot <- function(stanFit) {
  p2 <- mcmc_areas(
    as.array(stanFit), 
    pars = c("sigma_p", "sigma_pt","sigma_t","sigma_e"),
    prob = 0.8, # 80% intervals
    prob_outer = 0.95, 
    point_est = "median"
  )
  return(p2)
}


### CODE ###

# check for model_TCGApriors.rds, must add to prevent crash
checkFile <- "model_TCGApriors.rds"
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
emptyFit <- stan(file="model_TCGApriors.stan", data = stanDat1, chains = 0)

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
#   CpGscore <- log2(posterior[,1]/posterior[,2])
#   fitSumm <- rbind(fitSumm,rep(median(CpGscore),length(mInd)))
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
# CpGscore_C <- parData$'7'[,1]
# #rm(parData)
# gc()

# run in serial
# want to collect data of random effect coefficients (bcd), fixed effects (betaT,mu), and sigmas
# going to make it run the site and present the data within the same for loop, for data efficiency
ptm <- proc.time()
for(i in 1:nsites) {
  print(paste("site:", i,"/",nsites))
  
  data <- site(siteInds[i])
  stanFit <- stanfit_model(data)
  
  # take summary values for first 6 params (mu, nu, sigmas x4) and last param (lp__)
  fitSumm <- summary(stanFit)$summary[c(1:6,dim(summary(stanFit)$summary)[1]), c(1, 2, 6, 9, 10)]
  
  # compute posterior log-ratio: log2(sigmaP/sigmaT), take integral>0, mean, median
  posterior <- as.matrix(stanFit,pars=c("sigma_p","sigma_t"))
  logPTratio <- log2(posterior[,1]/posterior[,2])
  CpGscore <- c(mean(logPTratio > 0), mean(logPTratio), median(logPTratio))
  fitSumm <- rbind(fitSumm,CpGscore)
  
  # print summary of parameters (mu, nu, sigmas, lp) for this site
  print(summary(stanFit)$summary[c(1:6,dim(summary(stanFit)$summary)[1]),c(1, 2, 6, 9, 10)])
  
  # plot beta and fits at this site for all patients, using mean estimates for b and c coefficients
  p1 <- fitplot(data, stanFit)

  # plot posteriors for this site
  # variables: mu, nu, sigma_p, sigma_pt, sigma_t, sigma_e, lp__
  posterior <- as.array(stanFit)
  p2 <- mcmc_areas(
    posterior,
    pars = c("mu","nu"),
    prob = 0.8, # 80% intervals
    prob_outer = 0.95,
    point_est = "mean"
  ) + ggtitle(paste("Site",siteInds[i],": Fixed Effect Posteriors"))
  p3 <- mcmc_areas(
    posterior,
    pars = c("sigma_p","sigma_t","sigma_pt","sigma_e"),
    prob = 0.8, # 80% intervals
    prob_outer = 0.95,
    point_est = "mean"
  ) + ggtitle(paste("Site",siteInds[i],": Sigma Posteriors"))

  p4 <- traceplot(stanFit,pars=c("mu","nu"),inc_warmup=TRUE)
  p5 <- traceplot(stanFit,pars=c("sigma_e", "sigma_p", "sigma_pt","sigma_t"),inc_warmup=TRUE)

  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  
}
proc.time() - ptm

print(paste("Completed run on sites:", paste(siteInds, collapse=" ")))