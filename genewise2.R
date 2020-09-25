##### Gene-wise analysis
# UPDATED
# Kevin Murgas

library(tidyverse)
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")
library(here)

# load annotated data, Stan results
load(here("myFA_bulkonly.Rdata"))
load(here("FullResultsTCGA_relaxgamma3_postsum3.Rdata"))


med_vals <-
  tibble(
    Name = FullAnnotation[, 1],
    mu = mu_full$p50,
    betaT = betaT_full$p50,
    sigmaE = sigmaE_full$p50,
    sigmaP = sigmaP_full$p50,
    sigmaPT = sigmaPT_full$p50,
    sigmaT = sigmaT_full$p50,
    prob1 = prob1_full,
    prob2 = 1 - prob2_full
  ) %>% mutate(
    logP = log(sigmaP),
    logPT = log(sigmaPT),
    logT = log(sigmaT)
  )

tibble(FullAnnotation[,c(1,2,3,5,6,8)]) %>% rename(Name = IlmnID) %>%
  left_join(by="Name", med_vals)

# first take gene averages of every parameter (and # sites)




#### test code for plotting overlaid histograms ####

library(tidyverse)
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")
librray(here)

prepped_data <- read_csv(here("data-stan_rg3ps2", "prepped_data.csv"))
genemask = !is.na(prepped_data$UCSC_RefGene_Name)

fBG <- function(x, a, m1, s1, m2, s2) a*dnorm(x,m1,s1) + (1-a)*dnorm(x,m2,s2) # bimodal gaussian

post_prior_histogram <- function(postvals,mask,Rhat,priorfun,bw,xl,titlestr) {
  # stacked histogram colored by gene-associated sites
  # and prior curve overlaid
  plot_df <- data.frame(val=postvals, gene_assoc=mask)
  p1 <- ggplot(plot_df, aes(x=val)) +
    geom_histogram(aes(y=..count../length(mask),fill=gene_assoc), binwidth=bw) +
    stat_function(fun=function(x) bw*priorfun(x),geom="line",color="black") +
    ggtitle(titlestr) + xlim(xl) + theme_classic() + theme(legend.position = "none") +
    xlab(titlestr) + ylab("Freq.")
  
  # add overlay for Rhat > 1.1 (error sites to be filtered)
  errinds = which(Rhat > 1.1)
  err_df <- data.frame(val=postvals[errinds])
  p1 <- p1 + geom_density(data=err_df,mapping=aes(x=val,y=0.1*..density..),color="red")
  
  p1
}
post_prior_histogram(mu_full$p50, genemask, lp_full$Rhat,
                     function(x) fBG(x,0.3557,-3.0675, 0.75214,1.1290, 1.39197),
                     0.1, c(-5,5), "Mu")
post_prior_histogram(log(sigmaT_full$p50),genemask,lp_full$Rhat, 
                     function(x) 0.5*dgamma(exp(x), 3.3643, 10.430),
                     0.1, c(-5,5), "log(SigmaT)") + theme(legend.position = c(0.8, 0.6))



##### test for patient-specific analysis #####
# take top 50-100 genes or pathways
# use log transformed beta values for PWD (?)
# then plot a heat-map
# 3 forms:
#  1 = pairwise distance of each patients tumors
#  2 = same as above but normalized by PWD of all normal samples
#  3 = PWD of all tumors normalized by PWD of normal

library(gtools)

# Load raw methylation data
load(here("myFA_bulkonly.Rdata"))
permed_data <- read_csv(here("data-stan_rg3ps3", "permed_data.csv"))

# pull out all samples beta values
all_beta <- FullAnnotation[,10:57]

# extract patient and tissue labels
patientLabel <- substr(colnames(all_beta), 1, 1)
patientLabel[10:12] <- "K*"
tissueLabel <- substr(colnames(all_beta), 2, 2)
tissueLabel[!(tissueLabel %in% "N")] <- "T"   #replace A/B/D/M with T

pats = patientLabel %>% sort %>% unique
npat = length(pats)

permed_data <- read_csv(here("data-stan_rg3ps3", "permed_data.csv"))
#gened_data <- read_csv(here("data-stan_rg3ps3", "gened_data.csv"))

enrichedGenes <- permed_data %>%
  #filter(score_median_v3 > quantile(score_median_v3, 0.95, na.rm = TRUE)) %>%
  filter(score_median_v3 > 0.95) %>%
  select(UCSC_RefGene_Name, score_median, score_median_v3, n) %>%
  rename(Gene=UCSC_RefGene_Name, Pctile=score_median_v3)

ngene = 50 # using top 50 genes
ngene = nrow(enrichedGenes)
genes = enrichedGenes %>% top_n(-ngene,Pctile) %>% pull(Gene)
genes = genes[1:ngene]

#pwd <- function(x) {
#  mean(dist(x, method = "manhattan"))
#}
pwd <- function(x) {
  mean(dist(x))
}
pwd2 <- function(x) {
  abs(diff(x))
}

FullAnnotation_gened <- FullAnnotation %>%
  arrange(CHR, MAPINFO) %>%
  separate_rows(UCSC_RefGene_Name) %>% 
  distinct()
all_beta <- FullAnnotation_gened[,10:57]

signature_mat1 = matrix(,ngene,npat) # raw tumor pwd
for (igene in 1:ngene) {
  # each gene, get gene inds, calculate PWD for the tumor samples
  print(paste0("Gene: ",igene,"/",ngene))
  geneind = which(FullAnnotation_gened$UCSC_RefGene_Name %in% genes[igene])
  
  # calculate PWD of each tumor
  tum_pwd <- numeric(npat)
  for (ipat in 1:npat) {
    # each patient, get tumor sample inds
    patind <- which(patientLabel == pats[ipat] & tissueLabel == "T")
    tum_pwd[ipat] <- all_beta[geneind,patind] %>% apply(1,pwd2) %>% mean
  }
  
  # calculate PWD of all normal
  #norm_pwd <- all_beta[geneind,tissueLabel == "N"] %>% logit %>% apply(1, pwd) %>% mean
  
  signature_mat1[igene,] = tum_pwd
  #signature_mat2[igene,] = tum_pwd/norm_pwd
}

df1 = as.data.frame(signature_mat1)
rownames(df1) <- genes
colnames(df1) <- pats
write_csv(df1, here("data-stan_rg3ps3", "patmat1_nologit.csv"))

df1 <- read_csv(here("data-stan_rg3ps3", "patmat1_nologit.csv"))
df2 <- read_csv(here("data-stan_rg3ps3", "patmat2_nologit.csv"))
#read_csv(here("data-stan_rg3ps3", "patmat1_nologit.csv"))
#read_csv(here("data-stan_rg3ps3", "patmat2_nologit.csv"))

signature_mat1 %>% as.tibble %>% `colnames<-`(pats) %>% gather %>%
  mutate(gene=rep(genes,npat)) %>% rename(Patient=key,PWD=value) %>%
  ggplot + geom_tile(aes(x=gene,y=Patient,fill=PWD)) +
  theme(axis.text.x=element_text(angle=90,hjust=1))

df1 = as.data.frame(signature_mat1)
rownames(df1) <- genes
colnames(df1) <- pats
library(pheatmap)
pheatmap(
  #mat=signature_mat1,
  mat=t(df1),
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "test"
)
mat_cluster_cols <- hclust(dist(t(df1)))
plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")


library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dendsort)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

df1 <- read_csv(here("data-stan_rg3ps3", "patmat1_nologit.csv"))
rownames(df1) <- genes
mat_breaks <- quantile_breaks(as.matrix(df1), n = 11)
clust_pats <- sort_hclust(hclust(dist(t(df1))))
clust_genes <- sort_hclust(hclust(dist(df1)))
pheatmap(
  mat=df1,
  color = inferno(length(mat_breaks) - 1),
  breaks=mat_breaks,
  cluster_cols = clust_pats,
  cluster_rows = clust_genes,
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = "Top 1118 Genes: Tumor PWD"
)

df1 <- read_csv(here("data-stan_rg3ps3", "patmat1_logit.csv"))
rownames(df1) <- genes
clust_pats <- sort_hclust(hclust(dist(t(df1))))
clust_genes <- sort_hclust(hclust(dist(df1)))
pheatmap(
  mat=df1,
  cluster_cols = clust_pats,
  cluster_rows = clust_genes,
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = "Top 1118 Genes: Tumor PWD (logit)"
)

df1 <- read_csv(here("data-stan_rg3ps3", "patmat1_nologit.csv"))
rownames(df1) <- genes
clust_pats <- sort_hclust(hclust(dist(t(df1))))
clust_genes <- sort_hclust(hclust(dist(df1)))
pheatmap(
  mat=df1,
  cluster_cols = clust_pats,
  cluster_rows = clust_genes,
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = "Top 1118 Genes: Tumor PWD (no logit)"
)


# binary approach
# select each patients genes with lowest 25% PWD
df1 <- read_csv(here("data-stan_rg3ps3", "patmat1_nologit.csv"))
rownames(df1) <- genes
df_bin = as.data.frame(df1)
for (ipat in 1:ncol(df_bin)) {
  patdat = df_bin %>% pull(ipat)
  df_bin[,ipat] = patdat < quantile(patdat,0.25)
}
df_bin = 1*(df_bin)
clust_pats <- sort_hclust(hclust(dist(t(df_bin))))
clust_genes <- sort_hclust(hclust(dist(df_bin)))
pheatmap(
  mat=df_bin,
  cluster_cols = clust_pats,
  cluster_rows = clust_genes,
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = "Top 1118 Genes: Tumor PWD (no logit)"
)

# 3 blocks of genes reactome





##### test for GSAA #####
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")
library(here)
# Load raw methylation data
load(here("myFA_bulkonly.Rdata"))
permed_data <- read_csv(here("data-stan_rg3ps3", "permed_data.csv"))

# load in the C5 gene sets from MSigDB
library(GSA)
gmtfile <- file.path("/Users/kevinmurgas/Documents/Data+ project/EPIC data","c2.cp.v7.1.symbols.gmt")
genesets <- GSA.read.gmt(gmtfile)
nsets <- length(genesets[[1]])

# Procedure:
# 1. Order and rank all genes by evidence (currently using log-P/T-ratio integral > 0)
# 2. For each gene set:
#     Create two CDFs/running sums:
#     1. In-set genes: increment by the evidence magnitude or RANK
#     2. Out-of-set genes: increment by 1
#     add 0 to other sum if gene in/not in
#     Normalize by max at the end to get CDFs
#     Take modified K-S statistic, by measuring the difference in the max and min differences of the two CDFs
# 3. Each gene set now has a value to rank enrichment
#  KS statistic
# 4. For significance (FDR), need to compare to KS statistic
# if in set: X = sqrt((N-G)/G) ... with G = # genes in set, N = # genes total
# if not in set: XN=−sqrt(G/(N-G))

# make a data.frame for gene sets and scores
# columns for gene set name, score, num_genes, and skip based on filter <10 and >500
nsets <- length(genesets[[1]])
GSEA_df <- data.frame(GeneSet=genesets[[2]][1:nsets], score=numeric(nsets),
                      num_gene = unlist(lapply(genesets[[1]][1:nsets],length)),
                      num_gene_cut = numeric(nsets))

# have to order all genes first by rank of conservation (conserved genes come first)
genenames = permed_data$UCSC_RefGene_Name
ngene = length(genenames)

# reorder by evidence rank
ind = order(-permed_data$score_median)
#evidence = permed_data$score_median[ind]
evidence = seq(ngene,1,-1)
genenames = genenames[ind]

ptm <- proc.time()
for (iset in 1:nsets) {
  # take only genes that we have in full gene list
  set <- unlist(genesets[[1]][iset]) %>% intersect(genenames)
  GSEA_df$num_gene_cut[iset] = length(set)
  
  if (length(set)<10 | length(set)>500) {next} # skip if num_gene less than 10 or greater than 500
  cat(sprintf("\nRunning GSEA on set %i/%i: %s (%i genes)",iset,nsets,genesets[[2]][iset],length(set)))
  
  # compute running sum (as PDF first)
  runsum_in = numeric(ngene)
  runsum_out = numeric(ngene)
  for (i in 1:ngene) {
    if (genenames[i] %in% set) {
      # evidence = permed_data$median, or ranks
      runsum_in[i] = evidence[i]
      runsum_out[i] = 0
    } else {
      runsum_in[i] = 0
      runsum_out[i] = 1
    }
  }
  # convert to CDFs
  runsum_in <- cumsum(runsum_in)/sum(runsum_in)
  runsum_out <- cumsum(runsum_out)/sum(runsum_out)
  
  # compute K-S-like statistic
  v <- runsum_in - runsum_out
  score <- max(c(0,v)) + min(c(0,v))
  GSEA_df$score[iset] <- score
}
proc.time()-ptm

# make table of top scores
GSEA_df %>% arrange(-score) %>% select(GeneSet,score,num_gene) %>% slice(1:50) %>% gt()

# create bootstraps: 1000 random draws for each length of geneset
# first need to optimize
# ideas for optimizing:
# 1. remove for loop, use a faster string matching algorithm to create a binary in/out mask
# 2. use a single running sum variable, scale each by correct factor (for in and out)

# optimized version
genenames = permed_data$UCSC_RefGene_Name
ngene = length(genenames)
# reorder by evidence rank
genenames = genenames[order(-permed_data$score_median)]
evidence = seq(ngene,1,-1)
#evidence = abs(ngene/2 - evidence) # this is done in GSVA to upweight both tails...???
ptm <- proc.time()
for (iset in 1:nsets) {
  # take only genes that we have in full gene list
  set <- unlist(genesets[[1]][iset]) %>% intersect(genenames)
  GSEA_df$num_gene_cut[iset] = length(set)
  
  if (length(set)<10 | length(set)>500) {next} # skip if num_gene less than 10 or greater than 500
  #cat(sprintf("\nRunning GSEA on set %i/%i: %s (%i genes)",iset,nsets,genesets[[2]][iset],length(set)))
  
  # compute running sum (as PDF first)
  # create the random walk, with correct scaling for in/out genes
  #runsum = numeric(ngene)
  inmask = genenames %in% set
  #ranktot = sum(evidence[inmask]) # calculate sum of evidence in set for scaling
  runsum = evidence/sum(evidence[inmask]) * inmask - 1/(ngene-length(set)) * !inmask
  
  # convert to CDF, compute K-S-like statistic
  v <- cumsum(runsum)
  score <- max(c(0,v)) + min(c(0,v))
  GSEA_df$score[iset] <- score
}
proc.time()-ptm


# calculate bootstraps
setlengths <- unique(GSEA_df$num_gene_cut)
setlengths <- sort(setlengths[which(!(setlengths < 10 | setlengths > 500))])
# run thru each set length
# store 1000 bootstrapped scores using completely random gene sets
# later will map genesets to this bootstrap?

shuffle_genesets <- function(n_geneset, n_rep, data) {
  #replicate(n = n_rep, expr = sample(data, as.numeric(n_genes), TRUE))
  bootstrap = numeric(n_rep)
  for (ishuf in 1:n_rep) {
  # for random shuffles just make a random logical mask
  inds = sample(1:ngene,n_geneset, replace=F)
  inmask = logical(ngene)
  inmask[inds] = TRUE
  ranktot = sum(evidence[inmask]) # calculate sum of ranks in set ahead of time
  runsum = evidence/ranktot * inmask - 1/(ngene-length(set)) * !inmask
  
  # convert to CDF, compute K-S-like statistic
  v <- cumsum(runsum)
  score <- max(c(0,v)) + min(c(0,v))
  bootstrap[ishuf] <- score
  }
  bootstrap
}

names(setlengths) <- setlengths
ptm<-proc.time()
res <- lapply(setlengths, shuffle_genesets, n_rep = 1000, data = evidence)
proc.time()-ptm
plot(setlengths, lapply(res, mean) %>% unlist)
# map enrichment scores to get percentiles on bootstrap null distribution using upper tail probability
pvals <- map2_dbl(GSEA_df$score, GSEA_df$num_gene_cut, ~ mean(res[[as.character(.y)]] > .x))

# Benjamini-Hochberg step-up correction procedure
# this matches r p.adjust("BH")
ptemp <- pvals[!is.na(pvals)]
m <- length(ptemp)
BHp <- sort(ptemp,decreasing=T,na.last=TRUE)*m/seq(m,1)
BHp[order(ptemp,decreasing=T)] <- cummin(BHp)
q <- pvals*NA
q[!is.na(pvals)] <- BHp

GSEA_df$pval <- pvals
GSEA_df$FDR <- q
GSEA_df$padj <- p.adjust(pvals,method="BH")



##### differential methylation with minfi dmpfinder #####
library(tidyverse)
library(gtools)
library(ggplot2)
rename <- dplyr::rename

# get data
all_beta <- FullAnnotation[,10:57]
rownames(all_beta)<-FullAnnotation$IlmnID
tissueLabel <- substr(colnames(all_beta), 2, 2)
tissueLabel[!(tissueLabel %in% "N")] <- "T"

dmps_logit <- dmpFinder(logit(as.matrix(all_beta)), tissueLabel, type="categorical")
#dmps_logit <- dmps_logit[order(as.numeric(rownames(dmps_logit))),] # reorder back to original site order
dmps_logit <- dmps_logit %>% rownames_to_column %>% dplyr::slice(match(FullAnnotation$IlmnID,rownames(dmps_logit)))

cor(dmps_logit$intercept, mu_full$p50)

inds <- seq(1,866091,100)
inds <- inds[dmps_logit$qval[inds]<0.05]
plot(dmps_logit$intercept[inds], mu_full$p50[inds])

df <- cbind(dmps_logit,
            (betaT_full %>% select(p50) %>% dplyr::rename(betaT=p50)),
            (mu_full %>% select(p50) %>% dplyr::rename(mu=p50)),
            (sigmaPT_full %>% select(p50) %>% dplyr::rename(sigmaPT=p50)),
            (sigmaE_full %>% select(p50) %>% dplyr::rename(sigmaE=p50))
            ) %>%
  mutate(DMPsig=qval<0.05, betaTsig=abs(betaT)>quantile(abs(betaT),0.95), f2=abs(betaT)/sigmaE) %>%
  dplyr::slice(seq(1,866091,100))
ggplot(df) + geom_point(aes(x=f2,y=f,color=DMPsig), alpha=0.05) + ylim(0,25)




