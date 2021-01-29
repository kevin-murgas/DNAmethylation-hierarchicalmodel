# script for estimating sigma priors from TCGA patient data
# loads and combines patient data (with annotation)
# then algebraically calculates prior distribution estimates for each site
# to ultimately produce a full prior distribution then fit with a parameterized gamma
# currently test version

library(gtools)
library(ggplot2)

### INITIALIZE ###

#log.file <- "logTCGA.txt"
#writeLines(c(""), log.file)

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
#nsites <- 5
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
  nu <- line[2]
  linedf <- data.frame(tInd = c(0,1), est = c(mu, (mu+nu)))
  p1 <- p1 + geom_line(data=linedf, aes(x=tInd, y=est))
  
  return(p1)
}

calcSigmaAll <- function (siteData, plotFlag=FALSE) {
  siteBeta <- siteData$beta
  
  # sigmaP: take mean of normal samples
  mu <- mean(siteBeta[!tumorLabel])
  a <- siteBeta[!tumorLabel]
  sigmaP <- var(a)
  
  # sigmaPT: take mean of patient's tumor samples and normalize by patient normal
  normBeta <- c()
  for(pat in unique(patLabel)){
    if(length(intersect(which(pat==patLabel),which(!tumorLabel)))){ # only if pat has normal
      normIndex <- intersect(which(pat==patLabel),which(!tumorLabel)) # normal sample index
      tumIndex <- intersect(which(pat==patLabel),which(!!tumorLabel)) # tumor sample indices
      if(length(tumIndex)>1){ tum<-mean(siteBeta[tumIndex]) }else{ tum<-siteBeta[tumIndex] }
      normBeta <- c(normBeta, tum-siteBeta[normIndex])
    }
  }
  nu <- mean(normBeta)
  b <- normBeta
  sigmaPT <- var(b)
  
  # sigmaT: all samples var method, normalize all tumor samples by patient tumor mean
  normBetaT <- c()
  for(pat in unique(patLabel)){
    if(length(intersect(which(pat==patLabel),which(!!tumorLabel)))>1){ # only if pat has multiple tumor
      tumIndex <- intersect(which(pat==patLabel),which(!!tumorLabel)) # tumor sample indices
      mean <- mean(siteBeta[tumIndex])
      normBetaT <- c(normBetaT, siteBeta[tumIndex] - mean)
    }
  }
  g <- normBetaT
  sigmaT <- var(g)
  
  if (plotFlag == TRUE){
  ### PLOTS ###
  # plot sigmaP with normal distribution overlay
  data <- data.frame(beta=siteBeta[!tumorLabel])
  pP <- ggplot(data=data,aes(x=beta)) + geom_density()
  dataGen <- data.frame(beta=rnorm(10000,mean=mu,sd=sqrt(sigmaP)))
  pP <- pP + geom_density(data=dataGen, aes(x=beta), color="red") + xlim(-5,5)
  pP <- pP + labs(title="Mu+sigmaP: Normal samples",xlab="beta") + theme(plot.title = element_text(size = 8))
  #print(pP)
  
  # plot sigmaPT+nu with overlay
  data <- data.frame(beta=normBeta)
  pP2 <- ggplot(data=data,aes(x=beta)) + geom_density()
  dataGen <- data.frame(beta=rnorm(10000,mean=nu,sd=sqrt(sigmaPT)))
  pP2 <- pP2 + geom_density(data=dataGen, aes(x=beta), color="red") + xlim(-5,5)
  pP2 <- pP2 + labs(title="nu+sigmaPT: Tumor samples norm by Pat Normal",xlab="beta") + theme(plot.title = element_text(size = 7))
  #print(pP2)
  
  # plot sigmaT with overlay
  data <- data.frame(beta=normBetaT)
  pP3 <- ggplot(data=data,aes(x=beta)) + geom_density()
  dataGen <- data.frame(beta=rnorm(10000,mean=0,sd=sqrt(sigmaT)))
  pP3 <- pP3 + geom_density(data=dataGen, aes(x=beta), color="red") + xlim(-5,5)
  pP3 <- pP3 + labs(title="sigmaT: Tumor samples norm by Pat Tum mean",xlab="beta") + theme(plot.title = element_text(size = 7))
  #print(pP3)
  
  # fit plot
  # plot beta and fits at this site for all patients, using mean estimates for a and b coefficients
  p1 <- ggplot() + geom_point(data = siteData, aes(x=tInd, y=beta, colour=patient))
  pats <- as.integer(factor(data$patient))
  npats <- max(pats)
  linedf <- data.frame(tInd = c(0,1), est = c(mu, (mu+nu)))
  p1 <- p1 + geom_line(data=linedf, aes(x=tInd, y=est)) + guides(color="none")
  p1 <- p1 + labs(title="Site Data + Mean Fit") + theme(plot.title = element_text(size = 10))
  #print(p1)
  
  p2<-ggarrange(p1, pP, pP2, pP3, labels=c("A","B","C","D"), ncol=2,nrow = 2)
  print(p2)
  }
  
  return(c(mu,nu,sigmaP,sigmaPT,sigmaT))
}

# function to calculate sigmaP, uses only normal samples
# mu is mean, a is diff from mu, sigmaP is sample variance
calcSigmaP <- function (siteBeta) {
  mu <- mean(siteBeta[!tumorLabel])
  a <- siteBeta[!tumorLabel] - mu
  sigmaP <- var(a)
  
  return(c(mu,sigmaP))
}

# function to calculate sigmaPT, uses normal and tumor samples
# first take tumor mean (if multiple), then normalize by subtracting patient normal
# then mean is nu, b is diff from nu, sigmaPT is sample variance
calcSigmaPT <- function (siteBeta) {
  normBeta <- c()
  for(pat in unique(patLabel)){
    if(length(intersect(which(pat==patLabel),which(!tumorLabel)))){ # only if pat has normal
      normIndex <- intersect(which(pat==patLabel),which(!tumorLabel)) # normal sample index
      tumIndex <- intersect(which(pat==patLabel),which(!!tumorLabel)) # tumor sample indices
      if(length(tumIndex)>1){ tum<-mean(siteBeta[tumIndex]) }else{ tum<-siteBeta[tumIndex] }
      normBeta <- c(normBeta, tum-siteBeta[normIndex])
    }
  }
  nu <- mean(normBeta)
  b <- normBeta
  sigmaPT <- var(b)
  
  return(c(nu,sigmaPT))
}

# function to calculate sigmaT, uses patients with multiple tumor samples
# normalize by subtracting patient tumor mean
# then these values are g, sigmaT is sample variance
calcSigmaT <- function (siteBeta) {
  # all samples var method
  normBeta <- c()
  for(pat in unique(patLabel)){
    if(length(intersect(which(pat==patLabel),which(!!tumorLabel)))>1){ # only if pat has multiple tumor
      tumIndex <- intersect(which(pat==patLabel),which(!!tumorLabel)) # tumor sample indices
      mean <- mean(siteBeta[tumIndex])
      normBeta <- c(normBeta, siteBeta[tumIndex] - mean)
    }
  }
  g <- normBeta
  sigmaT <- var(g)
  
  # # mean of pat var method
  # popVars <- c()
  # for (pat in unique(patLabel)) {
  #   if (!!(length(intersect( which(pat == patLabel), which(!!tumorLabel))) - 1)) { # only if pat has multiple tumor
  #     tumIndex <- intersect(which(pat == patLabel), which(!!tumorLabel)) # tumor sample indices
  #     mean <- 0.5*(siteBeta[tumIndex[1]]-siteBeta[tumIndex[2]])^2
  #     popVars <- c(popVars, mean)
  #   }
  # }
  # sigmaT <- mean(popVars)
  
  return(sigmaT)
}


### CODE ###
# intialize numeric vectors for each site estimate
sigmaP_ests <- numeric(dim(concatData)[1])
sigmaPT_ests <- sigmaP_ests
sigmaT_ests <- sigmaP_ests
mu_ests <- sigmaP_ests
nu_ests <- sigmaP_ests

# get initial sample counts for a,b,g
a <- length(which(!tumorLabel))
b <- length(which(!tumorLabel))
g <- 0
for (pat in unique(patLabel)) {
  if(length(intersect(which(pat==patLabel),which(!!tumorLabel)))>1) { g <- g+length(intersect(which(pat==patLabel),which(!!tumorLabel))) }
}

print(c(a,b,g))

ptm<-proc.time()
for(i in siteInds){
  print(paste("Analyzing site",i,"."))
  siteData <- site(i)
  
  temp <- calcSigmaAll(siteData, plotFlag=FALSE)
  mu_ests[i] <- temp[1]
  nu_ests[i] <- temp[2]
  sigmaP_ests[i] <- temp[3]
  sigmaPT_ests[i] <- temp[4]
  sigmaT_ests[i] <- temp[5]
  print(paste("Site ",i," completed."))
  
  # plot raw data with mu and nu
  #p1 <- fitplot(site(i),c(mu_ests[i],nu_ests[i]))
  #print(p1)
  
  #invisible(readline(prompt="Press [enter] to continue"))
}
print(proc.time()-ptm)

save(sigmaP_ests,sigmaPT_ests,sigmaT_ests,mu_ests,nu_ests, file="estimates_allpats.Rdata")

########## FITTING DISTRIBUTIONS ##########
setwd("/Users/kevinmurgas/Documents/Data+ project/TCGA data")
load(file="estimates_allpats.Rdata")

# remove all zeros, and take square root of variance (stored variable) to get sigma/stdev
sigmaP_ests <- sqrt(sigmaP_ests[which(!!sigmaP_ests)])
sigmaPT_ests <- sqrt(sigmaPT_ests[which(!!sigmaPT_ests)])
sigmaT_ests <- sqrt(sigmaT_ests[which(!!sigmaT_ests)])
mu_ests <- mu_ests[which(!!mu_ests)]
nu_ests <-nu_ests[which(!!nu_ests)]

# examine histograms of estimates
hist(sigmaP_ests,breaks=seq(0,ceiling(max(sigmaP_ests)),0.02),xlim=c(0,2),main="sigmaP, all patients")
hist(sigmaPT_ests,breaks=seq(-3.2,ceiling(max(sigmaPT_ests)),0.02),xlim=c(0,2),main="sigmaPT, all patients")
hist(sigmaT_ests,breaks=seq(0,ceiling(max(sigmaT_ests)),0.02),xlim=c(0,2),main="sigmaT, all patients")
hist(mu_ests,breaks=100,xlim=c(-6,6),main="mu, all patients")
hist(nu_ests,breaks=100,xlim=c(-5,5),main="nu, all patients")

# now fit the distributions
library(MASS)

# fit mu with bimodal sum of gaussians
f_sumGauss <- function(x, a, mu1, sig1, mu2, sig2) ( a*dnorm(x,mean=mu1,sd=sig1) + (1-a)*dnorm(x,mean=mu2,sd=sig2))
mu_sumGauss <- fitdistr( mu_ests, f_sumGauss, list(a=0.5,mu1=-3.5,sig1=1,mu2=3.5,sig2=1), lower=c(0.001,-10,0.001,-10,0.001) , upper=c(1,10,10,10,10))

# fit nu with Cauchy distribution
nu_cauchy <- fitdistr( nu_ests, dcauchy, list(location=0,scale=1))

# fit sigmas with gamma
sigmaP_gamma <- fitdistr( sigmaP_ests, dgamma, list(shape=2, rate=2) )
sigmaPT_gamma <- fitdistr( sigmaPT_ests, dgamma, list(shape=2, rate=2) )
sigmaT_gamma <- fitdistr( sigmaT_ests, dgamma, list(shape=2, rate=2) )


# plot histograms again, with overlay of fits
# mu
f <- function(x, a, m1, s1, m2, s2) a*dnorm(x,m1,s1) + (1-a)*dnorm(x,m2,s2)
ggplot(data.frame(p=mu_ests), aes(x=p)) + geom_density(binwidth=0.01) + theme_bw() + ggtitle("Mu estimates and fit") + stat_function(fun=function(x) f(x,0.3557,-3.0675, 0.75214,1.1290, 1.39197),geom="line",color="red") + xlim(-6,6) + annotate("text",label="0.3557*N(-3.068, 0.752) +\n 0.6443*N(1.129, 1.392)",x=0,y=0.025,col="red")
# nu
hist(nu_ests,breaks=100,xlim=c(-6,6),main="nu",freq=FALSE)
curve( dcauchy(x,-0.04051, 0.19330),n=1001, add=TRUE, col="red")
text(3,1,"Cauchy(-0.04051, 0.19330)",col="red")
# sigmaP
hist(sigmaP_ests,breaks=seq(0,ceiling(max(sigmaP_ests)),0.01),xlim=c(0,1),ylim=c(0,11),main=expression(sigma["P"]),freq=FALSE)
curve( dgamma(x,5.1989, 16.6192),n=1001, add=TRUE, col="blue")
text(0.7,5,expression(Gamma*"("*alpha*"=5.1989,"~beta*"=16.6192)"),col="blue")
# sigmaPT
hist(sigmaPT_ests,breaks=seq(-2.2,ceiling(max(sigmaPT_ests)),0.02),xlim=c(0,2),ylim=c(0,4),main=expression(sigma["PT"]),freq=FALSE)
curve( dgamma(x,3.5539, 5.5324),n=1001, add=TRUE, col="blue")
text(1,2.5,expression(Gamma*"("*alpha*"=3.5539,"~beta*"=5.5324)"),col="blue")
# sigmaT
hist(sigmaT_ests,breaks=seq(0,ceiling(max(sigmaT_ests)),0.01),xlim=c(0,1),ylim=c(0,14),main=expression(sigma["T"]),freq=FALSE)
curve( dgamma(x,6.8383, 25.7560),n=1001, add=TRUE, col="blue")
text(0.7,5,expression(Gamma*"("*alpha*"=6.8383,"~beta*"=25.7560)"),col="blue")


# calculate relaxed priors (for gammas)
# increase variance by 2 but keep same mode
a1 = 5.1989; b1 = 16.6192 # sigmaP
a1 = 3.5539; b1 = 5.5324 # sigmaPT
a1 = 6.8383; b1 = 25.7560 # sigmaT
f1 <- function(a2, a1, b1) (a2-1)*b1/(a1-1) # mode
f2 <- function(a2, a1, b1) sqrt(a2*b1^2/(10*a1)) # 3x variance
a2range = c(0:1000000)*0.00001
out1 <- f1(a2range,a1,b1) - f2(a2range,a1,b1)
a2 = a2range[which(abs(out1)==min(abs(out1)))]
b2 = f1(a2,a1,b1)
plot(a2,b2,col='green',xlim=c(0,5),xlab='alpha2',ylab='beta2')
curve(f1(x,a1,b1),col='blue',add=TRUE)
curve(f2(x,a1,b1),col='red',add=TRUE)
curve(dgamma(x, a1, b1),n=1001,xlim=c(0,2),col='blue')
curve(dgamma(x, a2, b2),n=1001,xlim=c(0,2),add=TRUE,col='red')


# Examine priors coming out of stan
# define null model for stan, here with prior for sigmaP
library(rstan)
model = "parameters {real p;}
model {p ~ gamma(1.20572, 0.10636);}"
# model {target += log_mix(0.3557,normal_lpdf(p|-3.0675,0.56571) , normal_lpdf(p|1.1290, 1.9376));}"
fit = stan(model_code = model, pars = c('p'), control=list(adapt_delta=0.99, max_treedepth=10), iter = 1000, chains = 4, warmup = 500, verbose=FALSE)
posterior_P<-as.matrix(fit)

ggmcmc::ggs_density(ggmcmc::ggs(fit))+theme_minimal()

### for mu
# plot histogram, bimodal curve fit, and MCMC prior data
f <- function(x, a, m1, s1, m2, s2) a*dnorm(x,m1,s1) + (1-a)*dnorm(x,m2,s2)
plotData <- data.frame(p=mu_ests)
p1 <- ggplot(plotData, aes(x=p)) + geom_density(binwidth=0.01) + xlim(-6,6) + ggtitle("Mu estimates and fit")
#p1 <- p1 + stat_function(fun=function(x) f(x,0.444,-2.7,1.02,2,1.09),geom="line",color="red") + xlim(-6,6)
p2 <- mcmc_areas(posterior_P, pars = c("p"), prob = 0.8) + ggtitle('Mu: Bimodal Gaussian prior from Stan') + xlim(-6,6)
p3 <- ggmcmc::ggs_density(ggmcmc::ggs(fit))+theme_minimal()+xlim(-6,6) + ggtitle("Chains")
p4 <- ggarrange(p1, p2, p3, labels=c("A","B","C"), ncol=1,nrow = 3)
print(p4)

p1<-ggplot(plotData, aes(x=p))+stat_function(fun=function(x) f(x,1,1),geom="line",color="red") + xlim(0,10)
print(p1)

### for sigmas
# plot histogram (from model results), inv-gamma curve fit, and MCMC prior data
f <- function(x, rho, a) 1/(a*gamma(rho)) * (a / (x))^(rho+1) * exp( - a/(x) )
plotData <- data.frame(p=sigmaP_Cfull3$median) #Cfull3$mean
p1 <- ggplot(plotData, aes(x=p)) + geom_density(binwidth=0.01) + xlim(0,1) + ggtitle("SigmaPT estimates and fit")
p1 <- p1 + stat_function(fun=function(x) dgamma(x, shape=1.20572, scale=0.10636),geom="line",color="red")
print(p1)
p2 <- mcmc_areas(posterior_P, pars = c("p"), prob = 0.8) + ggtitle('Mu: Bimodal Gaussian prior from Stan') + xlim(0,6)
p3 <- ggmcmc::ggs_density(ggmcmc::ggs(fit))+theme_minimal()+xlim(0,6) + ggtitle("Chains")
p4 <- ggarrange(p1, p2, p3, labels=c("A","B","C"), ncol=1,nrow = 3)
print(p4)
