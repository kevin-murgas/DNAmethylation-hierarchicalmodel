# script for extracting variance from Stan complex model (model3)
# input taskID corresponding to 10000 site chunks starting at site (N-1)*10000+1
# Now incorporating parallel processing by RStan methods

library(doParallel)
library(rstan)
library(gtools)

rstan_options(auto_write = FALSE)
numCPUs = as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
print(paste("Cores: ", numCPUs))

# load data
load("myFA_bulkonly.Rdata")

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
  siteInds <- siteInds[1:6091]
}
nsites <- length(siteInds)

print(paste("Making chunk of nsites =", nsites))
chunk <- FullAnnotation[siteInds,]
rm(FullAnnotation)
gc()


### FXNS ###
# Function used to read in data from each site
# changed to now use a "chunk" instead of FullAnnotation, for speed and memory
# requires index 1-10000 instead of siteIndex
site <- function (ind) {
  # extract CpG site xx to start
  temp <- chunk[ind,]
  
  # here we use all patient samples, excluding glands
  indices <-  c(10:57) # if using bulk only data file
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
  
  return(work)
}

# Using Model with TCGA priors (previously model 3)
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

# Function for combining parallel runs of each fixef/sigma
comb <- function(x, ...) {
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}


### CODE ###

print("Generating data.frames")

mu_C <-
  data.frame(
    mean = numeric(nsites),
    sem = numeric(nsites),
    p50 = numeric(nsites),
    n_eff = numeric(nsites),
    Rhat = numeric(nsites)
  )
nu_C <- mu_C
sigmaP_C <- mu_C
sigmaPT_C <- mu_C
sigmaT_C <- mu_C
sigmaE_C <- mu_C
lp_C <- mu_C

# inds for extracting data (mean, sem, median, n_eff, Rhat)
mInd <- c(1, 2, 6, 9, 10)

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

# run in parallel via doParallel
cl <- makeCluster(numCPUs)
registerDoParallel(cl)
print(paste("Cores registered:",getDoParWorkers()))
print(paste("Backend type:",getDoParName()))
print("Starting foreach loop")
ptm <- proc.time()
parData <- foreach(i = iter(1:nsites), .combine = 'comb', .multicombine = TRUE, .packages=c("gtools","rstan")) %dopar% {
  if (i %% 1000 == 1) {
    sink(log.file, append = TRUE)
    print(paste("site:", i, "/", nsites))
    sink()
  }
  
  data <- site(i)
  stanFit <- stanfit_model(data)
  
  # take summary values for first 6 params (mu, nu, sigmas x4) and last param (lp__)
  fitSumm <- summary(stanFit)$summary[c(1:6,dim(summary(stanFit)$summary)[1]), mInd]
  
  # compute posterior log-ratio: log2(sigmaP/sigmaT), take integral>0, mean, median
  posterior <- as.matrix(stanFit,pars=c("sigma_p","sigma_t"))
  logPTratio <- log2(posterior[,1]/posterior[,2])
  CpGscore <- c(mean(logPTratio > 0), mean(logPTratio), median(logPTratio))
  fitSumm <- rbind(fitSumm,CpGscore)
  
  rm(data, stanFit, posterior, logPTratio, CpGscore)
  gc()
  
  split(fitSumm, rownames(fitSumm))
}
proc.time() - ptm
stopCluster(cl)
sink(log.file, append=TRUE)
print("All sites complete, saving data.")
sink()

mu_C[,] <- parData$mu
nu_C[,] <- parData$nu
sigmaP_C[,] <- parData$sigma_p
sigmaPT_C[,] <- parData$sigma_pt
sigmaT_C[,] <- parData$sigma_t
sigmaE_C[,] <- parData$sigma_e
lp_C[,] <- parData$lp__
CpGscore <- parData$CpGscore[,1:3]
rm(parData)
gc()

print(paste("Completed run, now saving"))
# save data
save(mu_C, nu_C, sigmaP_C, sigmaPT_C, sigmaT_C, sigmaE_C, lp_C, CpGscore, file = out.file)
