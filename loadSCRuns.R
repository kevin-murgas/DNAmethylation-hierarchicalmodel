# load StanCRun data
# produce full size data frames for mu,nu,sigmas (P,PT,T,E)
# read in individual save files containing 10000 site runs
# store 10000 site chunks in corresponding positions in full data frame
# repeat for each completed run

library(ggplot2)
library(reshape2)

### LOADING FOR STAN MODEL RESULTS

# Set working directory to normal data wd then folder with Results files from StanCParallel.R on all sites
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data/Stan Results Archive/FinalRuns/ResultsTCGA_gamma")

mu_full <-
  data.frame(
    mean = numeric(866091),
    sem = numeric(866091),
    p50 = numeric(866091),
    n_eff = numeric(866091),
    Rhat = numeric(866091)
  )

nu_full <- mu_full
sigmaP_full <- mu_full
sigmaPT_full <- mu_full
sigmaT_full <- mu_full
sigmaE_full <- mu_full
lp_full <- mu_full

CpGscore_full <- data.frame(prob = numeric(866091),
                         mean = numeric(866091),
                         median = numeric(866091))

nsites <- 10000
filesPresent <- (1:87)
fileNames <- dir(pattern = "StanCParResults_*")
allInds <- numeric()
for (i in filesPresent) {
  print(paste("Chunk:",i))
  load(fileNames[fileNames == paste("StanCParResults_",i,".Rdata", sep="")])
  inds <- (1:nsites) + (i - 1) * nsites
  if (i==87) {
    inds <- inds[1:6091]
  }
  allInds <- append(allInds, inds)
  mu_full[inds, ] <- mu_C
  nu_full[inds, ] <- nu_C
  sigmaP_full[inds, ] <- sigmaP_C
  sigmaPT_full[inds, ] <- sigmaPT_C
  sigmaT_full[inds, ] <- sigmaT_C
  sigmaE_full[inds, ] <- sigmaE_C
  lp_full[inds, ] <- lp_C
  CpGscore_full[inds, ] <- CpGscore
}
save(mu_full, nu_full, sigmaP_full, sigmaT_full, sigmaPT_full, sigmaE_full,
      lp_full, CpGscore_full, file = "FullResultsTCGA_gamma.Rdata")

