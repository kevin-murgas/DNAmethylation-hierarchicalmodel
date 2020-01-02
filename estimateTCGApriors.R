# script for estimating sigma priors from TCGA patient data
# loads and combines patient data (with annotation)
# then algebraically calculates prior distribution estimates for each site
# to ultimately produce a full prior distribution then fit with a parameterized gamma
# currently test version

library(gtools)
library(ggplot2)

### INITIALIZE ###

log.file <- "logTCGA.txt"
writeLines(c(""), log.file)

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/TCGA data")

# Yanlin's working directory
#setwd("D:/DataPlus2017/Data")

# Load and concatenate all available data
# need to: search all folders for data, load each one just the beta and label
fileList <- list.files(pattern = "_hg38.txt$", recursive = TRUE)
fileList <- fileList[-grep("tion27", fileList)] # filters out 27K files
# also check that A6-2672 is not utilized (only tumor B in 450K)

# intialize data frame
ptmTot <- proc.time()
print("Initializing data frame with annotations")
temp <- read.table(fileList[1],header=TRUE,sep="\t",fill=TRUE)
concatData <- as.data.frame.matrix(temp[,-2])

print("Now adding patient data")
for(f in fileList) {
  ptm <- proc.time()
  print(paste("File: ",which(f==fileList),"/",length(fileList)))
  dataLabel <- substr(f,regexpr("TCGA",f),regexpr("TCGA",f)+15)
  print(paste("Label: ", dataLabel))
  temp <- read.table(f,header=TRUE,sep="\t",fill=TRUE)
  concatData[dataLabel] <- temp[,2]
  print(proc.time() - ptm)
}
rm(temp)
gc()

print("Total time to load data:")
print(proc.time()-ptmTot)

save(concatData, file="testdata_allpats.Rdata")
load("testdata_allpats.Rdata") # start here to load concatenated patient data

# create labels
# start by making vectors for labels of patient and normal vs tumor
dataLabel <- colnames(concatData)[-1:-10]
patLabel <- substr(dataLabel,9,12)
tumorLabel <- 1-as.numeric(substr(dataLabel,14,14))
sideLabel <- substr(dataLabel,16,16)
sideLabel[grep("C",sideLabel)] <- "B" # simply re-labels any C tumor with B (should be only one case)
sideLabel[!tumorLabel] <- "N"

# find sites with NA
# need to check that all samples are not NA? or skip with any NA? or maybe at least half
# currently going to remove any site with NA
NAsites <- c()
for(i in dataLabel){
  NAsites <- union(NAsites,which(is.na(concatData[,i])))
}
print(paste("# sites with at least one NA: ",length(NAsites)))

# randomly select 5 sites in 1-485577 (excluding sites with any NA samples)
#nsites <- 10000
#siteInds <- sample(setdiff(1:485577,NAsites),nsites)

# alternatively, use the entire set
siteInds <- setdiff(1:485577,NAsites)
nsites <- length(siteInds)

### FXNS ###
# Extract single site data
site <- function (site_no) {
  # extract CpG site xx to start, exclude first 10 columns
  temp <- concatData[site_no,-1:-10]
  
  # already have dataLabel, patLabel, tumorLabel, sideLabel defined in initialization
  tissueLabel <- sideLabel
  tissueLabel[tissueLabel %in% c("A","B")] <- "T"   #replace A and B with T
  tumorIndicator <- 1*(tissueLabel=="T")
  
  work <- data.frame(beta=logit(t(temp)), patient = patLabel, tissue=tissueLabel, side=sideLabel, tInd=tumorIndicator)
  colnames(work)[1] <- "beta"
  return(work)
}

# Plot raw data and fits
fitplot <- function(data, line) {
  # plot beta and fits at this site for all patients, using mean estimates for b and c coefficients
  p1 <- ggplot() + geom_point(data = data, aes(x=tInd, y=beta, colour=patient))
  pats <- as.integer(factor(data$patient))
  npats <- max(pats)
  mu <- line[1]
  betaT <- line[2]
  linedf <- data.frame(tInd = c(0,1), est = c(mu, (mu+betaT)))
  p1 <- p1 + geom_line(data=linedf, aes(x=tInd, y=est))
  
  return(p1)
}

# function to calculate sigmaP, uses only normal samples
# mu is mean, b is diff from mu, sigmaP is sample variance
calcSigmaP <- function (siteBeta) {
  mu <- mean(siteBeta[!tumorLabel])
  b <- siteBeta[!tumorLabel] - mu
  sigmaP <- var(b)
  return(c(mu,sigmaP))
}

# function to calculate sigmaPT, uses normal and tumor samples
# first take tumor mean (if multiple), then normalize by subtracting patient normal
# then mean is betaT, c is diff from betaT, sigmaPT is sample variance
calcSigmaPT <- function (siteBeta) {
  normBeta <- c()
  for(pat in unique(patLabel)){
    if(length(intersect(which(pat==patLabel),which(!tumorLabel)))){ # only if pat has normal
      normIndex <- intersect(which(pat==patLabel),which(!tumorLabel)) # normal sample index
      tumIndex <- intersect(which(pat==patLabel),which(!!tumorLabel)) # tumor sample indices
      if(length(tumIndex)-1){ tum<-mean(siteBeta[tumIndex]) }else{ tum<-siteBeta[tumIndex] }
      normBeta <- c(normBeta, tum-siteBeta[normIndex])
    }
  }
  betaT <- mean(normBeta)
  c <- normBeta - betaT
  sigmaPT <- var(c)
  return(c(betaT,sigmaPT))
}

# function to calculate sigmaT, uses patients with multiple tumor samples
# #normalize by subtracting patient tumor mean
# #then these values are d, sigmaT is sample variance
# now uses mean of population variance for each patient with 2 samples
calcSigmaT <- function (siteBeta) {
  # # all samples var method
  # normBeta <- c()
  # for(pat in unique(patLabel)){
  #   if(!!(length(intersect(which(pat==patLabel),which(!!tumorLabel)))-1)){ # only if pat has multiple tumor
  #     tumIndex <- intersect(which(pat==patLabel),which(!!tumorLabel)) # tumor sample indices
  #     mean <- mean(siteBeta[tumIndex])
  #     normBeta <- c(normBeta, siteBeta[tumIndex] - mean)
  #   }
  # }
  # d <- normBeta
  # sigmaT <- var(d)
  
  # mean of pat var method
  popVars <- c()
  for (pat in unique(patLabel)) {
    if (!!(length(intersect( which(pat == patLabel), which(!!tumorLabel))) - 1)) { # only if pat has multiple tumor
      tumIndex <- intersect(which(pat == patLabel), which(!!tumorLabel)) # tumor sample indices
      mean <- 0.5*(siteBeta[tumIndex[1]]-siteBeta[tumIndex[2]])^2
      popVars <- c(popVars, mean)
    }
  }
  sigmaT <- mean(popVars)
  
  return(sigmaT)
}


### CODE ###
# intialize numeric vectors for each site estimate
sigmaP_ests <- numeric(dim(concatData)[1])
sigmaPT_ests <- sigmaP_ests
sigmaT_ests <- sigmaP_ests
mu_ests <- sigmaP_ests
betaT_ests <- sigmaP_ests

# get initial sample counts for b,c,d
b <- length(which(!tumorLabel))
c <- length(which(!tumorLabel))
d <- 0
for (pat in unique(patLabel)) {
  if(!!(length(intersect(which(pat==patLabel),which(!!tumorLabel)))-1)) { d <- d+length(intersect(which(pat==patLabel),which(!!tumorLabel))) }
}

print(c(b,c,d))

ptm<-proc.time()
for(i in siteInds){
  print(paste("Analyzing site",i,"."))
  siteData <- site(i)
  
  temp <- calcSigmaP(siteData$beta)
  mu_ests[i] <- temp[1]
  sigmaP_ests[i] <- temp[2]
  temp <- calcSigmaPT(siteData$beta)
  betaT_ests[i] <- temp[1]
  sigmaPT_ests[i] <- temp[2]
  sigmaT_ests[i] <- calcSigmaT(siteData$beta)
  sigmaPT_ests[i] <- sigmaPT_ests[i] - sigmaT_ests[i]  # subtract out sigmaT from sigmaPT
  print(paste("Site ",i," completed."))
  
  # plot raw data with mu and beta
  #p1 <- fitplot(site(i),c(mu_ests[i],betaT_ests[i]))
  #print(p1)
  
  #invisible(readline(prompt="Press [enter] to continue"))
}
print(proc.time()-ptm)

save(sigmaP_ests,sigmaPT_ests,sigmaT_ests,mu_ests,betaT_ests, file="estimates_allpats.Rdata")
load(file="estimates_allpats.Rdata")

logtemp <- log(sigmaP_ests[which(!!sigmaP_ests)])
skewnessLog <- 3*(mean(logtemp) - median(logtemp))/sd(logtemp)

# remove all zeros
sigmaP_ests <- sigmaP_ests[which(!!sigmaP_ests)]
sigmaPT_ests <- sigmaPT_ests[which(!!sigmaPT_ests)]
sigmaT_ests <- sigmaT_ests[which(!!sigmaT_ests)]
mu_ests <- mu_ests[which(!!mu_ests)]
betaT_ests <- betaT_ests[which(!!betaT_ests)]

# examine histograms of estimates
hist(sigmaP_ests,breaks=seq(0,ceiling(max(sigmaP_ests)),0.02),xlim=c(0,2),main="sigmaP, all patients")
hist(sigmaPT_ests,breaks=seq(-3.2,ceiling(max(sigmaPT_ests)),0.02),xlim=c(-1.5,3),main="sigmaPT, all patients")
hist(sigmaT_ests,breaks=seq(0,ceiling(max(sigmaT_ests)),0.02),xlim=c(0,2),main="sigmaT, all patients")
hist(mu_ests,breaks=100,xlim=c(-6,6),main="mu, all patients")
hist(betaT_ests,breaks=100,xlim=c(-5,5),main="betaT, all patients")

# repeat for log transform of sigmas (examine skew to determine distribution)
# calculate skew first
skewP <- 3*(mean(log(sigmaP_ests)) - median(log(sigmaP_ests)))/ sd(log(sigmaP_ests))
skewPT <- 3*(mean(log(sigmaPT_ests)) - median(log(sigmaPT_ests)))/ sd(log(sigmaPT_ests))
skewT <- 3*(mean(log(sigmaT_ests)) - median(log(sigmaT_ests)))/ sd(log(sigmaT_ests))
hist(log(sigmaP_ests),breaks=seq(-10,10,0.2),xlim=c(-10,10),main="log(sigmaP), all patients") + text(x=4,y=20000, paste("Skew =", format(skewP, nsmall=3)))
hist(log(sigmaPT_ests),breaks=seq(-10,10,0.2),xlim=c(-10,10),main="log(sigmaPT), all patients") + text(x=4,y=20000, paste("Skew =", format(skewPT, nsmall=3)))
hist(log(sigmaT_ests),breaks=seq(-10,10,0.2),xlim=c(-10,10),main="log(sigmaT), all patients") + text(x=4,y=20000, paste("Skew =", format(skewT, nsmall=3)))


# look at sigmaP/sigmaT
PTratio <- sigmaP_ests/sigmaT_ests
hist(PTratio,breaks=seq(0,7000,0.2),xlim=c(0,10),main="PTratio, all patients")
hist(log(PTratio),breaks=seq(-10,10,0.2),xlim=c(-10,10),main="log(PTratio), all patients")
skewLPT <- 3*(mean(log(PTratio)) - median(log(PTratio)))/ sd(log(PTratio))

# now fit the distributions
library(MASS)
# fit sigmas with inv-gamma
f <- function(x, rho, a) 1/(a*gamma(rho)) * (a / (x))^(rho+1) * exp( - a/(x) )
sigmaP_invgamma <- fitdistr( sigmaP_ests, f, list(rho=1, a=0.1) )
sigmaPT_invgamma <- fitdistr( sigmaPT_ests, f, list(rho=1, a=0.1) )
sigmaT_invgamma <- fitdistr( sigmaT_ests, f, list(rho=1, a=0.1) )

# mu will be large normal (K=10^4)
# instead of fitting, use 1/var = precision for betaT
tau <- 1/var(betaT_ests)


# plot histograms again, with overlay of fits
hist(sigmaP_ests,breaks=seq(0,ceiling(max(sigmaP_ests)),0.01),xlim=c(0,1),ylim=c(0,11),main=expression(sigma["P"]*","~Gamma^-1*"("*alpha*"=1.7969,"~beta*"=0.10645)"),freq=FALSE)
curve( f(x,sigmaP_invgamma$estimate[1],sigmaP_invgamma$estimate[2]),n=1001, add=TRUE, col="red")

hist(sigmaPT_ests,breaks=seq(0,ceiling(max(sigmaPT_ests)),0.02),xlim=c(0,2),ylim=c(0,4),main=expression(sigma["PT"]*","~Gamma^-1*"("*alpha*"=0.90825,"~beta*"=0.14678)"),freq=FALSE)
curve( f(x,sigmaPT_invgamma$estimate[1],sigmaPT_invgamma$estimate[2]),n=1001, add=TRUE, col="red")

hist(sigmaT_ests,breaks=seq(0,ceiling(max(sigmaT_ests)),0.01),xlim=c(0,1),ylim=c(0,14),main=expression(sigma["T"]*","~Gamma^-1*"("*alpha*"=1.7004,"~beta*"=0.07478)"),freq=FALSE)
curve( f(x,sigmaT_invgamma$estimate[1],sigmaT_invgamma$estimate[2]),n=1001, add=TRUE, col="red")

hist(betaT_ests,breaks=100,xlim=c(-5,5),main="betaT, norm(0,sd=0.62121) tau=2.5914",freq=FALSE)
curve( dnorm(x,mean=0,sd=1/sqrt(tau)), add=TRUE, col="red")


# more test code
plot(density(sigmaP_ests),col="blue",xlim=c(0,1),ylim=c(0,11),main=expression(sigma["P"]*","~Gamma^-1))
hist(sigmaP_ests,breaks=seq(0,ceiling(max(sigmaP_ests)),0.01),add=TRUE,freq=FALSE)
curve( f(x,sigmaP_invgamma$estimate[1],sigmaP_invgamma$estimate[2]),n=1001, add=TRUE, col="red")

plot(density(sigmaPT_ests),col="blue",xlim=c(0,2),ylim=c(0,4),main=expression(sigma["PT"]*","~Gamma^-1))
hist(sigmaPT_ests,breaks=seq(0,ceiling(max(sigmaPT_ests)),0.02),add=TRUE,freq=FALSE)
curve( f(x,sigmaPT_invgamma$estimate[1],sigmaPT_invgamma$estimate[2]),n=1001, add=TRUE, col="red")

plot(density(sigmaT_ests),col="blue",xlim=c(0,1),ylim=c(0,14),main=expression(sigma["T"]*","~Gamma^-1))
hist(sigmaT_ests,breaks=seq(0,ceiling(max(sigmaT_ests)),0.01),add=TRUE,freq=FALSE)
curve( f(x,sigmaT_invgamma$estimate[1],sigmaT_invgamma$estimate[2]),n=1001, add=TRUE, col="red")
