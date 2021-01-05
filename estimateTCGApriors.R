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
  betaT <- line[2]
  linedf <- data.frame(tInd = c(0,1), est = c(mu, (mu+betaT)))
  p1 <- p1 + geom_line(data=linedf, aes(x=tInd, y=est))
  
  return(p1)
}

calcSigmaAll <- function (siteData, plotFlag=FALSE) {
  siteBeta <- siteData$beta
  
  # sigmaP: take mean of normal samples
  mu <- mean(siteBeta[!tumorLabel])
  b <- siteBeta[!tumorLabel]
  sigmaP <- var(b)
  
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
  betaT <- mean(normBeta)
  c <- normBeta
  sigmaPT <- var(c)
  
  # sigmaT: all samples var method, normalize all tumor samples by patient tumor mean
  normBetaT <- c()
  for(pat in unique(patLabel)){
    if(length(intersect(which(pat==patLabel),which(!!tumorLabel)))>1){ # only if pat has multiple tumor
      tumIndex <- intersect(which(pat==patLabel),which(!!tumorLabel)) # tumor sample indices
      mean <- mean(siteBeta[tumIndex])
      normBetaT <- c(normBetaT, siteBeta[tumIndex] - mean)
    }
  }
  d <- normBetaT
  sigmaT <- var(d)
  
  if (plotFlag == TRUE){
  ### PLOTS ###
  # plot sigmaP with normal distribution overlay
  data <- data.frame(beta=siteBeta[!tumorLabel])
  pP <- ggplot(data=data,aes(x=beta)) + geom_density()
  dataGen <- data.frame(beta=rnorm(10000,mean=mu,sd=sqrt(sigmaP)))
  pP <- pP + geom_density(data=dataGen, aes(x=beta), color="red") + xlim(-5,5)
  pP <- pP + labs(title="Mu+sigmaP: Normal samples",xlab="beta") + theme(plot.title = element_text(size = 8))
  #print(pP)
  
  # plot sigmaPT+betaT with overlay
  data <- data.frame(beta=normBeta)
  pP2 <- ggplot(data=data,aes(x=beta)) + geom_density()
  dataGen <- data.frame(beta=rnorm(10000,mean=betaT,sd=sqrt(sigmaPT)))
  pP2 <- pP2 + geom_density(data=dataGen, aes(x=beta), color="red") + xlim(-5,5)
  pP2 <- pP2 + labs(title="betaT+sigmaPT: Tumor samples norm by Pat Normal",xlab="beta") + theme(plot.title = element_text(size = 7))
  #print(pP2)
  
  # plot sigmaT with overlay
  data <- data.frame(beta=normBetaT)
  pP3 <- ggplot(data=data,aes(x=beta)) + geom_density()
  dataGen <- data.frame(beta=rnorm(10000,mean=0,sd=sqrt(sigmaT)))
  pP3 <- pP3 + geom_density(data=dataGen, aes(x=beta), color="red") + xlim(-5,5)
  pP3 <- pP3 + labs(title="sigmaT: Tumor samples norm by Pat Tum mean",xlab="beta") + theme(plot.title = element_text(size = 7))
  #print(pP3)
  
  # fit plot
  # plot beta and fits at this site for all patients, using mean estimates for b and c coefficients
  p1 <- ggplot() + geom_point(data = siteData, aes(x=tInd, y=beta, colour=patient))
  pats <- as.integer(factor(data$patient))
  npats <- max(pats)
  linedf <- data.frame(tInd = c(0,1), est = c(mu, (mu+betaT)))
  p1 <- p1 + geom_line(data=linedf, aes(x=tInd, y=est)) + guides(color="none")
  p1 <- p1 + labs(title="Site Data + Mean Fit") + theme(plot.title = element_text(size = 10))
  #print(p1)
  
  p2<-ggarrange(p1, pP, pP2, pP3, labels=c("A","B","C","D"), ncol=2,nrow = 2)
  print(p2)
  }
  
  return(c(mu,betaT,sigmaP,sigmaPT,sigmaT))
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
      if(length(tumIndex)>1){ tum<-mean(siteBeta[tumIndex]) }else{ tum<-siteBeta[tumIndex] }
      normBeta <- c(normBeta, tum-siteBeta[normIndex])
    }
  }
  betaT <- mean(normBeta)
  c <- normBeta
  sigmaPT <- var(c)
  
  return(c(betaT,sigmaPT))
}

# function to calculate sigmaT, uses patients with multiple tumor samples
# normalize by subtracting patient tumor mean
# then these values are d, sigmaT is sample variance
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
  d <- normBeta
  sigmaT <- var(d)
  
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
betaT_ests <- sigmaP_ests

# get initial sample counts for b,c,d
b <- length(which(!tumorLabel))
c <- length(which(!tumorLabel))
d <- 0
for (pat in unique(patLabel)) {
  if(length(intersect(which(pat==patLabel),which(!!tumorLabel)))>1) { d <- d+length(intersect(which(pat==patLabel),which(!!tumorLabel))) }
}

print(c(b,c,d))

ptm<-proc.time()
for(i in siteInds){
  print(paste("Analyzing site",i,"."))
  siteData <- site(i)
  
  # temp <- calcSigmaP(siteData$beta)
  # mu_ests[i] <- temp[1]
  # sigmaP_ests[i] <- temp[2]
  # temp <- calcSigmaPT(siteData$beta)
  # betaT_ests[i] <- temp[1]
  # sigmaPT_ests[i] <- temp[2]
  # sigmaT_ests[i] <- calcSigmaT(siteData$beta)
  # sigmaPT_ests[i] <- sigmaPT_ests[i] - sigmaT_ests[i]  # subtract out sigmaT from sigmaPT
  
  temp <- calcSigmaAll(siteData, plotFlag=FALSE)
  mu_ests[i] <- temp[1]
  betaT_ests[i] <- temp[2]
  sigmaP_ests[i] <- temp[3]
  sigmaPT_ests[i] <- temp[4]
  sigmaT_ests[i] <- temp[5]
  # sigmaPT_ests[i] <- sigmaPT_ests[i] - sigmaT_ests[i]  # subtract out sigmaT from sigmaPT?
  print(paste("Site ",i," completed."))
  
  # plot raw data with mu and beta
  #p1 <- fitplot(site(i),c(mu_ests[i],betaT_ests[i]))
  #print(p1)
  
  #invisible(readline(prompt="Press [enter] to continue"))
}
print(proc.time()-ptm)

save(sigmaP_ests,sigmaPT_ests,sigmaT_ests,mu_ests,betaT_ests, file="estimates_allpats.Rdata")

########## FITTING DISTRIBUTIONS ##########
setwd("/Users/kevinmurgas/Documents/Data+ project/TCGA data")
load(file="estimates_allpats.Rdata")

# remove all zeros, and take square root of variance (stored variable) to get sigma/stdev
sigmaP_ests <- sqrt(sigmaP_ests[which(!!sigmaP_ests)])
sigmaPT_ests <- sqrt(sigmaPT_ests[which(!!sigmaPT_ests)])
sigmaT_ests <- sqrt(sigmaT_ests[which(!!sigmaT_ests)])
mu_ests <- mu_ests[which(!!mu_ests)]
betaT_ests <- betaT_ests[which(!!betaT_ests)]

# examine histograms of estimates
hist(sigmaP_ests,breaks=seq(0,ceiling(max(sigmaP_ests)),0.02),xlim=c(0,2),main="sigmaP, all patients")
hist(sigmaPT_ests,breaks=seq(-3.2,ceiling(max(sigmaPT_ests)),0.02),xlim=c(0,2),main="sigmaPT, all patients")
hist(sigmaT_ests,breaks=seq(0,ceiling(max(sigmaT_ests)),0.02),xlim=c(0,2),main="sigmaT, all patients")
hist(mu_ests,breaks=100,xlim=c(-6,6),main="mu, all patients")
hist(betaT_ests,breaks=100,xlim=c(-5,5),main="betaT, all patients")

# now fit the distributions
library(MASS)

# fit mu with bimodal sum of gaussians
f_sumGauss <- function(x, a1, mu1, sig1, a2, mu2, sig2) ( a1/(sig1*sqrt(2*pi)) * exp(-0.5*((x-mu1)/sig1)^2) ) + ( a2/(sig2*sqrt(2*pi)) * exp(-0.5*((x-mu2)/sig2)^2) )
f_sumGauss2 <- function(x, a, mu1, sig1, mu2, sig2) ( a*dnorm(x,mean=mu1,sd=sig1) + (1-a)*dnorm(x,mean=mu2,sd=sig2))
f_sumGauss3 <- function(x, mu1, sig1, mu2, sig2) ( 0.5*dnorm(x,mean=mu1,sd=sig1) + 0.5*dnorm(x,mean=mu2,sd=sig2))
mu_sumGauss <- fitdistr( mu_ests, f_sumGauss2, list(a1=0.5,mu1=-3.5,sig1=1,a2=0.5,mu2=3.5,sig2=1), lower=c(0.001,-10,0.001,0.001,-10,0.001) , upper=c(1,10,10,1,10,10))
mu_sumGauss <- fitdistr( mu_ests, f_sumGauss3, list(mu1=3,sig1=1,mu2=-3,sig2=1) )

# alternatively fit mu with logit transformed beta distributions

# fit betaT with Cauchy distribution
betaT_cauchy <- fitdistr( betaT_ests, dnorm, list(mean=0,sd=1))

# fit sigmas with inv-gamma (rho=alpha, a=beta)
f <- function(x, rho, a) 1/(a*gamma(rho)) * (a / (x))^(rho+1) * exp( - a/(x) )
sigmaP_invgamma <- fitdistr( sigmaP_ests, f, list(rho=1, a=0.1) )
sigmaPT_invgamma <- fitdistr( sigmaPT_ests, f, list(rho=1, a=0.1) )
sigmaT_invgamma <- fitdistr( sigmaT_ests, f, list(rho=1, a=0.1) )
# also try gamma?
sigmaP_gamma <- fitdistr( sigmaP_ests, dgamma, list(shape=2, rate=2) )
sigmaPT_gamma <- fitdistr( sigmaPT_ests, dgamma, list(shape=2, rate=2) )
sigmaT_gamma <- fitdistr( sigmaT_ests, dgamma, list(shape=2, rate=2) )


# plot histograms again, with overlay of fits
# mu
f <- function(x, a, m1, s1, m2, s2) a*dnorm(x,m1,s1) + (1-a)*dnorm(x,m2,s2)
ggplot(data.frame(p=mu_ests), aes(x=p)) + geom_density(binwidth=0.01) + theme_bw() + ggtitle("Mu estimates and fit") + stat_function(fun=function(x) f(x,0.3557,-3.0675, 0.75214,1.1290, 1.39197),geom="line",color="red") + xlim(-6,6) + annotate("text",label="0.3557*N(-3.068, 0.752) +\n 0.6443*N(1.129, 1.392)",x=0,y=0.025,col="red")
hist(mu_ests,breaks=100,xlim=c(-6,6),main="mu estimates",freq=FALSE)
curve( f_sumGauss3(x,mu_sumGauss$estimate[1],mu_sumGauss$estimate[2],mu_sumGauss$estimate[3],mu_sumGauss$estimate[4]), n=1001, add=TRUE, col="red")
# betaT
hist(betaT_ests,breaks=100,xlim=c(-6,6),main="betaT",freq=FALSE)
curve( dcauchy(x,-0.04051, 0.19330),n=1001, add=TRUE, col="red")
curve( dnorm(x,mean=-0.057128,sd=0.621209), add=TRUE, col="blue")
text(3,1,"Cauchy(-0.04051, 0.19330)",col="red")
text(3,0.8,"Normal(-0.057128, 0.621209)",col="blue")
# sigmaP
f <- function(x, rho, a) 1/(a*gamma(rho)) * (a / (x))^(rho+1) * exp( - a/(x) )
hist(sigmaP_ests,breaks=seq(0,ceiling(max(sigmaP_ests)),0.01),xlim=c(0,1),ylim=c(0,11),main=expression(sigma["P"]),freq=FALSE)
curve( f(x,6.3641, 1.6631),n=1001, add=TRUE, col="red")
curve( dgamma(x,5.1989, 16.6192),n=1001, add=TRUE, col="blue")
#curve( f(x,sigmaP_invgamma$estimate[1],sigmaP_invgamma$estimate[2]),n=1001, add=TRUE, col="red")
#curve( dgamma(x,sigmaP_gamma$estimate[1],sigmaP_gamma$estimate[2]),n=1001, add=TRUE, col="blue")
text(0.7,6,expression(Gamma^-1*"("*alpha*"=6.3641,"~beta*"=1.6631)"),col="red")
text(0.7,5,expression(Gamma*"("*alpha*"=5.1989,"~beta*"=16.6192)"),col="blue")
# sigmaPT
hist(sigmaPT_ests,breaks=seq(-2.2,ceiling(max(sigmaPT_ests)),0.02),xlim=c(0,2),ylim=c(0,4),main=expression(sigma["PT"]),freq=FALSE)
curve( f(x,3.2624, 1.5397),n=1001, add=TRUE, col="red")
curve( dgamma(x,3.5539, 5.5324),n=1001, add=TRUE, col="blue")
text(1,3,expression(Gamma^-1*"("*alpha*"=3.2624,"~beta*"=1.5397)"),col="red")
text(1,2.5,expression(Gamma*"("*alpha*"=3.5539,"~beta*"=5.5324)"),col="blue")
# sigmaT
hist(sigmaT_ests,breaks=seq(0,ceiling(max(sigmaT_ests)),0.01),xlim=c(0,1),ylim=c(0,14),main=expression(sigma["T"]),freq=FALSE)
curve( f(x,6.5392, 1.4894),n=1001, add=TRUE, col="red")
curve( dgamma(x,6.8383, 25.7560),n=1001, add=TRUE, col="blue")
text(0.7,6,expression(Gamma^-1*"("*alpha*"=6.5392,"~beta*"=1.4894)"),col="red")
text(0.7,5,expression(Gamma*"("*alpha*"=6.8383,"~beta*"=25.7560)"),col="blue")



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


# calculate relaxed priors (for gammas)
# increase variance by 2 but keep same mode
a1 = 6.8383; b1 = 25.7560
f1 <- function(a2, a1, b1) (a2-1)*b1/(a1-1) # mode
f2 <- function(a2, a1, b1) sqrt(a2*b1^2/(3*a1)) # 3x variance
a2range = c(0:1000000)*0.00001
out1 <- f1(a2range,a1,b1) - f2(a2range,a1,b1)
a2 = a2range[which(abs(out1)==min(abs(out1)))]
b2 = f1(a2,a1,b1)
plot(a2,b2,col='green',xlim=c(0,5),xlab='alpha2',ylab='beta2')
curve(f1(x,a1,b1),col='blue',add=TRUE)
curve(f2(x,a1,b1),col='red',add=TRUE)
curve(dgamma(x, a1, b1),n=1001,xlim=c(0,2),col='blue')
curve(dgamma(x, a2, b2),n=1001,xlim=c(0,2),add=TRUE,col='red')

# calculate relaxed priors (for inverse gammas)
# increase variance by 2 but keep same mode
dIG <- function(x, a, b) 1/(b*gamma(a))*(b/x)^(a+1)*exp(-b/x)
a1 = 2.5; b1 = 2.5
f1 <- function(a2, a1, b1) b1/(a1+1)*(a2+1) # mode
f2 <- function(a2, a1, b1) sqrt((100*b1^2*(a2-1)^2*(a2-2))/((a1-1)^2*(a1-2))) # 10x variance
a2range = c(1:10000000)*0.000001+2
out1 <- f1(a2range,a1,b1) - f2(a2range,a1,b1)
a2 = a2range[which(abs(out1)==min(abs(out1)))]
b2 = f1(a2,a1,b1)
plot(a2,b2,col='green',xlim=c(0,5),xlab='alpha2',ylab='beta2')
curve(f1(x,a1,b1),col='blue',add=TRUE)
curve(f2(x,a1,b1),col='red',add=TRUE)
curve(dIG(x, a1, b1),n=1001,xlim=c(0,2),col='blue')
curve(dIG(x, a2, b2),n=1001,xlim=c(0,2),add=TRUE,col='red')
curve(dIG(exp(x), a1, b1),n=1001,xlim=c(-5,2),col='blue')
curve(dIG(exp(x), a2, b2),n=1001,xlim=c(-5,2),add=TRUE,col='red')



# fit bimodal mu with EM algorithm
library(dplyr)
library(rearrange2)
e_step <- function(x, mu.vector, sd.vector, alpha.vector) {
  comp1.prod <- dnorm(x, mu.vector[1], sd.vector[1]) * alpha.vector[1]
  comp2.prod <- dnorm(x, mu.vector[2], sd.vector[2]) * alpha.vector[2]
  sum.of.comps <- comp1.prod + comp2.prod
  comp1.post <- comp1.prod / sum.of.comps
  comp2.post <- comp2.prod / sum.of.comps
  
  sum.of.comps.ln <- log(sum.of.comps, base = exp(1))
  sum.of.comps.ln.sum <- sum(sum.of.comps.ln)
  
  list("loglik" = sum.of.comps.ln.sum,
       "posterior.df" = cbind(comp1.post, comp2.post))
}
m_step <- function(x, posterior.df) {
  comp1.n <- sum(posterior.df[, 1])
  comp2.n <- sum(posterior.df[, 2])
  
  comp1.mu <- 1/comp1.n * sum(posterior.df[, 1] * x)
  comp2.mu <- 1/comp2.n * sum(posterior.df[, 2] * x)
  
  comp1.var <- sum(posterior.df[, 1] * (x - comp1.mu)^2) * 1/comp1.n
  comp2.var <- sum(posterior.df[, 2] * (x - comp2.mu)^2) * 1/comp2.n
  
  comp1.alpha <- comp1.n / length(x)
  comp2.alpha <- comp2.n / length(x)
  
  list("mu" = c(comp1.mu, comp2.mu),
       "var" = c(comp1.var, comp2.var),
       "alpha" = c(comp1.alpha, comp2.alpha))
}
wait <- mu_ests
wait.kmeans <- kmeans(wait, 2)
wait.kmeans.cluster <- wait.kmeans$cluster
wait.df <- data.frame(x = wait, cluster = wait.kmeans.cluster)
ggplot(data=wait.df, aes(y=wait, x=x, color = factor(cluster))) + geom_point() +ylab("Values") +ylab("Data Point Number") +scale_color_discrete(name = "Cluster") +ggtitle("K-means Clustering")
wait.df %>%
  mutate(num = row_number()) %>%
  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  geom_point() +
  ylab("Values") +
  ylab("Data Point Number") +
  scale_color_discrete(name = "Cluster") +
  ggtitle("K-means Clustering")
wait.summary.df <- wait.df %>%
  group_by(cluster) %>%
  summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
wait.summary.df %>%
  select(cluster, mu, variance, std)
wait.summary.df <- wait.summary.df %>%
  mutate(alpha = size / sum(size))
wait.summary.df %>%
  select(cluster, size, alpha)
for (i in 1:100) {
  print(i)
  if (i == 1) {
    # Initialization
    e.step <- e_step(wait, wait.summary.df[["mu"]], wait.summary.df[["std"]],
                     wait.summary.df[["alpha"]])
    m.step <- m_step(wait, e.step[["posterior.df"]])
    cur.loglik <- e.step[["loglik"]]
    loglik.vector <- e.step[["loglik"]]
  } else {
    # Repeat E and M steps till convergence
    e.step <- e_step(wait, m.step[["mu"]], sqrt(m.step[["var"]]), 
                     m.step[["alpha"]])
    m.step <- m_step(wait, e.step[["posterior.df"]])
    loglik.vector <- c(loglik.vector, e.step[["loglik"]])
    
    loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
    if(loglik.diff < 1e-6) {
      break
    } else {
      cur.loglik <- e.step[["loglik"]]
    }
  }
}
loglik.vector
m.step
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}
data.frame(x = wait) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 0.1, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(m.step$mu[1], sqrt(m.step$var[1]), 
                            lam = m.step$alpha[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(m.step$mu[2], sqrt(m.step$var[2]), 
                            lam = m.step$alpha[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density") +
  xlab("Values") +
  ggtitle("Final GMM Fit")


# try with optim
library(MASS)
loglik <- function(mu,x) {
  sum(-mu[5]*dbeta(x,mu[1],mu[3],log=TRUE) - mu[6]*dbeta(x,mu[2],mu[4],log=TRUE))
}
out <- optim(par = c(2,8,8,2,1,1), fn=loglik,x=mu_ests,method = "L-BFGS-B",lower=c(0,0,0,0,0,0))

