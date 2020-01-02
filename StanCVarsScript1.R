# script for extracting variance from Stan model on random 5K sites
# part 1


library(rstan)
library(gtools)

# load data
load("myFA.Rdata")

args <- commandArgs(trailingOnly = TRUE)
N <- as.numeric(args[1])
out.file <- paste("Results_", args[1], ".Rdata", sep="")

# Function used to read in data from each site

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
  
  return(work)
}



# Function to build stan model for each site

stanfit3 <- function (dataset) {
  
  stanDat <- list(pID = as.integer(factor(dataset$patient)),
                  tInd = dataset$tInd,
                  N = nrow(dataset),
                  P = nlevels(dataset$patient),
                  y = dataset[,1])
  
  
  # Using Model 3: add intra-tumoral variances
  stanFit3 <- stan(file="model3.stan", data=stanDat, control = list(adapt_delta = 0.999))
  
  return(stanFit3=stanFit3)
}



# Function to get mode

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


# Choose 1k sites by N

siteInds <- ((N*1000-999):(N*1000))
if (N>866) { siteInds <- ((N*1000-999):866836) }
nsites <- length(siteInds)



betaT_C <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                      p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                      p97.5=numeric(nsites))

mu_C <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                   p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                   p97.5=numeric(nsites))

sigmaE_C <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                       p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                       p97.5=numeric(nsites))


sigmaP_C <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                       p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                       p97.5=numeric(nsites))


sigmaT_C <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                       p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                       p97.5=numeric(nsites))


sigmaPT_C <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                        p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                        p97.5=numeric(nsites))

count <- 1

for (i in indice_1) {
  data <- site(i)
  stan <- stanfit3(data)
  posterior <- as.array(stan)
  
  betaT_C[count,1] <- summary(stan)$summary[71,1]
  betaT_C[count,2] <- getmode(posterior[1:1000,1:4,71])
  betaT_C[count,3:7] <- summary(stan)$summary[71,4:8]
  
  mu_C[count,1] <- summary(stan)$summary[72,1]
  mu_C[count,2] <- getmode(posterior[1:1000,1:4,72])
  mu_C[count,3:7] <- summary(stan)$summary[72,4:8]
  
  sigmaE_C[count,1] <- summary(stan)$summary[73,1]
  sigmaE_C[count,2] <- getmode(posterior[1:1000,1:4,73])
  sigmaE_C[count,3:7] <- summary(stan)$summary[73,4:8]
  
  
  sigmaP_C[count,1] <- summary(stan)$summary[74,1]
  sigmaP_C[count,2] <- getmode(posterior[1:1000,1:4,74])
  sigmaP_C[count,3:7] <- summary(stan)$summary[74,4:8]
  
  sigmaPT_C[count,1] <- summary(stan)$summary[75,1]
  sigmaPT_C[count,2] <- getmode(posterior[1:1000,1:4,75])
  sigmaPT_C[count,3:7] <- summary(stan)$summary[75,4:8]
  
  sigmaT_C[count,1] <- summary(stan)$summary[76,1]
  sigmaT_C[count,2] <- getmode(posterior[1:1000,1:4,76])
  sigmaT_C[count,3:7] <- summary(stan)$summary[76,4:8]
  
  count <- count + 1
}


# save data
save(mu_C,betaT_C,sigmaP_C,sigmaT_C,
     sigmaPT_C,sigmaE_C,file=out.file)
