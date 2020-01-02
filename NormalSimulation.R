#################################################################################
######################### Normal Tissue Data Simulation #########################
#################################################################################


library(rstan)
library(gtools)
library(ggplot2)
library(bayesplot)
library(MASS)

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# Yanlin's working directory
#setwd("D:/DataPlus2017/Data")

# load data
load("myFA.Rdata")


# Extract single site data
site <- function (site_no) {
  # extract CpG site xx to start
  temp <- FullAnnotation[site_no,]
  
  # here we use all patient samples, excluding glands
  indices <- c(9:13, 15:39,46,47,57,58,72:75)
  patientLabel <- substr(colnames(temp[indices]),1,1)
  patientLabel[10:12] <- "K*"
  sideLabel <- substr(colnames(temp[indices]),2,2)
  tissueLabel <- sideLabel
  tissueLabel[tissueLabel %in% c("A","B")] <- "T"   #replace A and B with T
  tumorIndicator <- 1*(tissueLabel=="T")
  
  work <- data.frame(beta=logit(t(temp[indices])), patient = patientLabel, tissue=tissueLabel, side=sideLabel, tInd=tumorIndicator)
  colnames(work)[1] <- "beta"
  return(work)
}

# Using Model 3: add intra-tumoral variance

stanfit3 <- function (dataset) {
  
  stanDat <- list(pID = as.integer(factor(dataset$patient)),
                  tInd = dataset$tInd,
                  N = nrow(dataset),
                  P = nlevels(dataset$patient),
                  y = dataset[,1])
  
  
  
  stanFit3 <- stan(file="model3.stan", data=stanDat, control = list(adapt_delta = 0.999))
  
  return(stanFit3=stanFit3)
}

# Plot posterior distribution
posterior <- function(dataset) {
  stan <- stanfit3(dataset)
  mcmc_areas(
    as.array(stan), 
    pars = c("sigma_p", "sigma_pt","sigma_t","sigma_e"),
    prob = 0.8, # 80% intervals
    prob_outer = 0.95, 
    point_est = "median"
  )
}

# Randomly choose 5 sites
set.seed(555)
siteInds <- sample(1:866836,5)
# [1] 321382 765044 600495 465787 559918
nsites <- length(siteInds)


# Original posterior
data321383 <- site(siteInds[1])
data765044 <- site(siteInds[2])
data600495 <- site(siteInds[3])
data465787 <- site(siteInds[4])
data559918 <- site(siteInds[5])

posterior(data321383)
posterior(data765044)
posterior(data600495)
posterior(data465787)
posterior(data559918)


#### Simulate normal tissue data and redo the bayesian analysis
simulation <- function (dataset) {
  pat_tissue_unique <- unique(dataset[,c("patient","tissue")])
  dup <- (duplicated(unique(dataset[,c("patient","tissue")])[,"patient"]))
  pat_normal <- dataset$patient[((dataset$patient) %in% (pat_tissue_unique$patient[dup]))]
  pat_all <- unique(dataset$patient)
  fit <- fitdistr(dataset$beta[dataset$tissue=="N"],"normal")

  for (i in pat_all) {
    if (i %in% pat_normal == FALSE) {
      dataset <- rbind(dataset,
                          c(rnorm(1,fit$estimate[1],fit$estimate[2]),i,"N","N",0))
    }
  }

  dataset[,1] <- as.numeric(dataset[,1])
  dataset[,5] <- as.numeric(dataset[,5])

  return(dataset)
}

data321383 <- simulation(data321383)
data765044 <- simulation(data765044)
data600495 <- simulation(data600495)
data465787 <- simulation(data465787)
data559918 <- simulation(data559918)


posterior(data321383)
posterior(data765044)
posterior(data600495)
posterior(data465787)
posterior(data559918)
