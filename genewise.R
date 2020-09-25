##### Gene-wise analyses
# Kevin Murgas
# Looking at CpG sites within each individual gene
# Using full-scale (850K site) Stan model results

# Some analyses below (will) include:
# Gene Average of Sigma Parameters (PCA, clustering)
# Distance-based analysis (sigmas vs CpG distance in gene)
# F vs R strand methylation correlation
# Differential methylation validation with mu+betaT

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# load annotation, Stan results
load("myFA_bulkonly.Rdata")
load("FullResultsTCGA_relaxgamma3_postsum3.Rdata")

# load GeneInfo (gene names and regions for each site)
load("GeneInfo.Rdata")

### if GeneInfo not made yet ###

# make siteGenes, a list of genes at each site
fsplit <- function(str) unique(unlist(strsplit(str,split=";")))
siteGenes <- lapply(FullAnnotation$UCSC_RefGene_Name, fsplit) # takes ~38 sec

# make siteRegions, a label for regulatory region at each site
# use 5'>3' priority: TSS1500>TSS200>5UTR>1stExon>Body>ExonBnd>3UTR>blank
fsplitReg <- function(str) {
    worklist <- unique(unlist(strsplit(str, split = ";")))
    worklist[worklist == "3'UTR"] <- "3UTR"
    worklist[worklist == "5'UTR"] <- "5UTR"
    worklist[worklist == "5URT"] <- "5UTR"
    #if (length(worklist)==0) { worklist <- "blank"}
    return(unique(worklist))
  }
siteRegions <- lapply(FullAnnotation$UCSC_RefGene_Group, fsplitReg)

# pre-index all genes
geneInd <- function (geneStr) {
  strx = paste0("(?<=^|;)", geneStr, "(?=$|;)") # match string - look for start/end of string or semicolon on either side of match
  matchinds = grep(strx,FullAnnotation$UCSC_RefGene_Name,perl=TRUE); # check UCSC name
  return(matchinds)
}
ptm<-proc.time()
gene2sites <- lapply(uniqueGenes, geneInd) # runs in ~50 mins
proc.time()-ptm

# save these data as GeneInfo.Rdata
save(siteGenes, siteRegions, gene2sites, file = "GeneInfo.Rdata")
load("GeneInfo.Rdata")

# for each gene, find gene start site from Ensembl annotation database
# store the gene range, and direction, and use direction to find start site (or do that later)
library(EnsDb.Hsapiens.v75) # ensembl annotation with gene position (hg19 version)
edb <- EnsDb.Hsapiens.v75


# get unique genes/regions, number of genes per site, sites per gene
uniqueGenes <- unique(unlist(siteGenes))
ngenes <- length(uniqueGenes)
uniqueRegions <- unique(unlist(siteRegions))
nregions <- length(uniqueRegions)
lenSG <- unlist(lapply(siteGenes, length))
lenSR <- unlist(lapply(siteRegions, length))
lenGS <- unlist(lapply(gene2sites,length))

# make a dataframe of the 3 sigma vectors (take median of posterior)
sigma_df <- data.frame(sigmaP=sigmaP_full$p50,sigmaPT=sigmaPT_full$p50,sigmaT=sigmaT_full$p50)
logsig_df <- log(sigma_df)

##### Functions ######
# find indices corresponding to geneString in the list geneNames, now using regex for speed
geneInd <- function (geneStr) {
  strx = paste("(?<=^|;)", geneStr, "(?=$|;)",sep="") # match string - look for start/end of string or semicolon on either side of match
  matchinds = grep(strx,FullAnnotation$UCSC_RefGene_Name,perl=TRUE); # check UCSC name
  return(matchinds)
}
# version 2 has the RefGene names as a string array, not sure if faster
geneStrings <- FullAnnotation$UCSC_RefGene_Name
geneInd_regexp <- function (geneStr) {
  strx = paste("(?<=^|;)", geneStr, "(?=$|;)",sep="") # match string - look for start/end of string or semicolon on either side of match
  matchinds = grep(strx,geneStrings,perl=TRUE); # check UCSC name
  return(matchinds)
}



##### Make a dummy Methylset to get genomic locations (not working) #####




##### PCA on gene averages #####


##### CpG-distance analysis #####
# this will need to be parallelized and batch processed...
# goal is to loop thru every gene, make a vector of CpG distances (same length as sigma vectors)
# also store number of sites in a gene
#library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
#library(ExpressionAtlas)
#library(org.Hs.eg.db) # organism level (mostly gene aliases)
#library(GO.db) # gene ontology
library(EnsDb.Hsapiens.v75) # ensembl annotation with gene position (hg19 version)
edb <- EnsDb.Hsapiens.v75
columns(edb) # list columns of db
keytypes(edb) # list keys of db
geneKeys <- keys(edb, keytype="GENENAME")
geneKeys <- intersect(uniqueGenes,geneKeys)
edbData<-select(edb, keys=geneKeys, keytype="GENENAME", columns=c("GENEBIOTYPE", "GENESEQSTART", "GENESEQEND", "SEQNAME", "SEQSTRAND"))
edbData<- edbData[which(a$SEQNAME %in% c(as.character(1:22),"X","Y")),] # keep only rows which map to chromosomes not patches

ngenes <- length(uniqueGenes)
diffPos <- numeric(ngenes)
flagGene <- logical(ngenes)
for (i in 1:ngenes) {
  geneStr <- uniqueGenes[i]
  inds <- which(edbData$GENENAME==geneStr)
  print(paste("Gene:",i,"/",ngenes,":",geneStr,"(nind =",length(inds),all(edbData$SEQSTRAND[inds]==edbData$SEQSTRAND[inds[1]]),")"))
  #if (length(inds) == 0) { print("0 ind")}
  #if (length(inds) > 1) { print(paste("multiple inds, (",length(inds),") same strand?",all(edbData$SEQSTRAND[inds]==edbData$SEQSTRAND[inds[1]])))}
  if (length(inds) > 1 && !all(edbData$SEQSTRAND[inds]==edbData$SEQSTRAND[inds[1]])) {flagGene[i] <- TRUE}
  meanEdb <- mean(c(mean(edbData$GENESEQSTART[inds]),mean(edbData$GENESEQEND[inds])))
  meanCpG <- mean(FullAnnotation$MAPINFO[gene2sites[[which(uniqueGenes==geneStr)]]])
  #print(paste("Mean Pos (EDB):",meanEdb))
  #print(paste("Mean Pos (CpG):",meanCpG))
  diffPos[i] <- meanEdb-meanCpG
}
# bad: SNORA75=333, PIAS3=1290
# good multiple: TSEN34, SHANK2, LILRA6, MROH1, NFKBIL1, ATP6V1G2, DDX5, CDSN, PSORS1C1
# big diff: RALGAPA1P


# loop thru every gene to do: F vs R, distance
CpGdist <- numeric(dim(FullAnnotation)[1]) + NA
ngenes <- length(uniqueGenes)
Fmeth <- numeric(ngenes) + NA
Rmeth <- Fmeth
n<-numeric(ngenes)
allInds <- c()
ptm<-proc.time()
for (i in 1:ngenes) { # this loop takes about 50 min to run
  # find all sites in gene
  print(paste("Gene:",i,"/",ngenes))
  gInd <- geneInd(uniqueGenes[i])
  allInds <- c(allInds,gInd)
  # check that all sites on same chromosome and same strand dir
  # ? sites can be on F or R strands, not sure how to tell which is 5' for the gene....
  # just go with lowest MAPINFO as left-most for now....
  #n[i] <- length(gInd)
  #pos <- FullAnnotation$MAPINFO[gInd]
  #left <- min(pos)
  #dist <- pos - left
  #CpGdist[gInd] <- dist
  Fmeth[i] <- mean(mu_full$p50[gInd[FullAnnotation$Strand[gInd]=="F"]])
  Rmeth[i] <- mean(mu_full$p50[gInd[FullAnnotation$Strand[gInd]=="R"]])
}
print(proc.time()-ptm)

# save CpG distance analysis results?

# plot sigmas based on CpGdist bins
# every 100 bp up to 20000
binBounds = seq(0, 20000, 100)
nBin = length(binBounds) - 1
binCents = (binBounds[1:nBin] + binBounds[2:(nBin+1)])/2 #seq(50, 19950, 100)
binMeans = numeric(length(binCents))
binSDs = numeric(length(binCents))
for (i in c(1:length(binCents))) {
  inds = which((CpGdist > binBounds[i]) & (CpGdist < binBounds[i+1]))
  inds = intersect(inds, which(sigmaP_full$n_eff>100))
  binMeans[i] = mean(sigmaP_full$mean[inds])
  binSDs[i] = sd(sigmaP_full$mean[inds])
}
ggplot(data=data.frame(cent=binCents,mean=binMeans,sd=binSDs),aes(x=cent,y=mean,err=sd)) + geom_pointrange(ymin=binMeans-binSDs,ymax=binMeans+binSDs)+ xlab("Bin Center") + ylab("Bin Mean") + ggtitle("CpG-Distance Analysis")


# new gene score
# ( sigmaPT^2 + sigmaT^2 )/ sigmaP^2
gsc <- (sigmaPT_full$p50^2 + sigmaT_full$p50^2)/ (2*sigmaP_full$p50^2)
hist(gsc)
hist(log(gsc))
hist(log10(gsc))
sum(gsc<1)
sum(gsc<1)/length(gsc)








##### fixing up the distance analysis #####

# load
library(tidyverse)
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")
library(here)
#load(here("myFA_bulkonly.Rdata"))
#load(here("FullResultsTCGA_relaxgamma3_postsum3.Rdata"))
gened_data <- read_csv(here("data-stan_rg3ps3", "gened_data.csv"), col_types = list("CHR" = col_character()))
#gened_data <- read_csv(here("data-stan", "gened_data.csv"), col_types = list("CHR" = col_character()))

# load ensembl db
library(EnsDb.Hsapiens.v75) # ensembl annotation with gene position (hg19 version)
edb <- EnsDb.Hsapiens.v75
geneKeys <- keys(edb, keytype="GENENAME") %>% intersect(unique(gened_data$UCSC_RefGene_Name))
edbData<- AnnotationDbi::select(edb, keys=geneKeys, keytype="GENENAME",
                                columns=c("GENEBIOTYPE", "GENESEQSTART", "GENESEQEND", "SEQNAME", "SEQSTRAND")) %>%
  filter(SEQNAME %in% c(as.character(1:22),"X","Y")) %>%
  mutate(GENESTART = GENESEQSTART*(SEQSTRAND==1) + GENESEQEND*(SEQSTRAND==-1)) %>%
  select(-c(GENESEQSTART,GENESEQEND))

# calculate CpG distance based on Ensembl gene start position and direction
distanced_data <- gened_data %>%
  filter(!is.na(UCSC_RefGene_Name)) %>%
  filter(UCSC_RefGene_Name %in% geneKeys) %>%
  mutate(EDBCHR=edbData$SEQNAME[match(UCSC_RefGene_Name,edbData$GENENAME)],
         EDBSTART=edbData$GENESTART[match(UCSC_RefGene_Name,edbData$GENENAME)],
         EDBDIR=edbData$SEQSTRAND[match(UCSC_RefGene_Name,edbData$GENENAME)],
         EDBTYPE=edbData$GENEBIOTYPE[match(UCSC_RefGene_Name,edbData$GENENAME)]
         ) %>%
  filter(CHR==EDBCHR) %>%
  mutate(cpgdist = if_else(EDBDIR==1, MAPINFO-EDBSTART, EDBSTART-MAPINFO))

distance_length <- distanced_data %>%
  filter(!is.na(UCSC_RefGene_Name)) %>%
  #slice(grep("Body",distanced_data$UCSC_RefGene_Group)) %>%
  #filter(EDBTYPE=="protein_coding") %>%
  select(CHR, UCSC_RefGene_Name, cpgdist, c(sigmaE:sigmaT,score_median:score_reg)) %>%
  #select(CHR, UCSC_RefGene_Name, cpgdist, c(sigmaP:sigmaT)) %>%
  #rename(score_median=sigmaP) %>% mutate(score_median2=score_median, score_reg=score_median) %>%
  #mutate(score_median=sigmaP,score_median2=sigmaP, score_reg=sigmaP) %>%
  mutate(location_bin = cut(cpgdist, c(seq(-1500, 20000, by = 50)), include.lowest = TRUE)) %>%
  group_by(location_bin) %>%
  summarise(mean_sigmaE = mean(sigmaE), mean_sigmaP = mean(sigmaP), mean_sigmaPT = mean(sigmaPT), mean_sigmaT = mean(sigmaT),  mean_score_median = mean(score_median), mean_score_median2 = mean(score_median2), mean_score_reg = mean(score_reg), n_obs = n(), .groups='drop') %>%
  mutate(location_bin = as.numeric(str_extract(location_bin, "(?<=,)(.*)(?=])")))


# # Figure 1 (original method, no direction correction)
# distance_length <- gened_data %>%
#   filter(!is.na(UCSC_RefGene_Name)) %>%
#   select(CHR, MAPINFO, UCSC_RefGene_Name, c(sigmaP:sigmaT,score_median:score_reg)) %>%
#   group_by(UCSC_RefGene_Name, CHR) %>%
#   mutate(location = MAPINFO - MAPINFO[1]) %>%
#   mutate(location_bin = cut(location, c(-1, seq(0, 20000, by = 100)), include.lowest = TRUE)) %>%
#   group_by(location_bin) %>%
#   summarise(mean_sigmaP = mean(sigmaP), mean_sigmaPT = mean(sigmaPT), mean_sigmaT = mean(sigmaT),  mean_score_median = mean(score_median), mean_score_median2 = mean(score_median2), mean_score_reg = mean(score_reg), n_obs = n()) %>%
#   mutate(location_bin = as.numeric(str_extract(location_bin, "(?<=,)(.*)(?=])")))

distance_length %>%
  drop_na() %>%
  pivot_longer(mean_sigmaE:mean_sigmaT) %>% 
  #mutate(name = factor(name, c("mean_sigmaP","mean_sigmaPT","mean_sigmaT"), c("SigmaP","SigmaPT","SigmaT"))) %>%
  mutate(name = factor(name, c("mean_sigmaE","mean_sigmaP","mean_sigmaPT","mean_sigmaT"), c("SigmaE","SigmaP","SigmaPT","SigmaT"))) %>%
  sample_frac() %>%
  ggplot(aes(location_bin, value, color = name, fill = name)) +
  geom_point(key_glyph = draw_key_rect) +
  ggforce::facet_zoom(location_bin < 5000) +
  theme_light() +
  theme(legend.position = "bottom") +
  labs(title = NULL,
       fill = NULL,
       color = NULL,
       x = "Distance from gene start site (binned to 100 bp)",
       y = "Average Sigma Value", 
       caption = "Figure 1: Average Sigma Values for single CpGs as a function of position relative to gene start site.") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2")

# again for conservation scores
distance_length %>%
  drop_na() %>%
  pivot_longer(mean_score_median) %>%
  mutate(name = factor(name, c("mean_score_median"), c("Score_logPTratio_median"))) %>%
  sample_frac() %>%
  ggplot(aes(location_bin, value, color = name, fill = name)) +
  geom_point(key_glyph = draw_key_rect) +
  ggforce::facet_zoom(location_bin < 5000) +
  theme_light() +
  theme(legend.position = "bottom") +
  labs(title = NULL,
       fill = NULL,
       color = NULL,
       x = "Distance from gene start site (binned to 100 bp)",
       y = "Average Conservation Score", 
       caption = "Figure 1: Average Conservation Scores for single CpGs as a function of position relative to gene start site.") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2")

# and look at number of CpG in each bin
distance_length %>%
  drop_na() %>%
  pivot_longer(n_obs) %>%
  mutate(name = factor(name, c("n_obs"), c("n_obs"))) %>%
  sample_frac() %>%
  ggplot(aes(location_bin, value, color = name, fill = name)) +
  geom_point(key_glyph = draw_key_rect) +
  ggforce::facet_zoom(location_bin < 5000) +
  theme_light() +
  theme(legend.position = "bottom") +
  labs(title = NULL,
       fill = NULL,
       color = NULL,
       x = "Distance from gene start site (binned to 100 bp)",
       y = "Number of CpGs", 
       caption = "Figure 1: Frequency of single CpGs as a function of position relative to gene start site.") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2")
