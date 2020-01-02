################### Plots ############################

library(ggplot2)
library(dunn.test)

# here set the working directory that points to the data folder
# e.g. the folder with annotated data saved as "myFA.Rdata"
# please just comment out the directories already here, so each
# user can uncomment the setwd for them as we move code around

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# Yanlin's working directory
setwd("D:/DataPlus2017/Data")


### LOAD DATA
load("StanCFullResults.RData")

sigma <- data.frame(sigmaP_Cfull$p50,sigmaPT_Cfull$p50,sigmaT_Cfull$p50,sigmaE_Cfull$p50)


### Estimates Distribution
ggplot(mu_Cfull,aes(x=mu_Cfull$p50)) + geom_histogram(binwidth = 0.1) + labs(title="Distribution of mu",x="mu (median)")

ggplot(betaT_Cfull,aes(x=betaT_Cfull$p50)) + geom_histogram(binwidth = 0.1) + labs(title="Distribution of betaT",x="betaT (median)")

ggplot(sigmaP_Cfull,aes(x=sigmaP_Cfull$p50)) + geom_histogram(binwidth = 0.05) + labs(title="Distribution of sigmaP",x="sigmaP (median)") + coord_cartesian(ylim=c(0,175000),xlim=c(0,3))

ggplot(sigmaPT_Cfull,aes(x=sigmaPT_Cfull$p50)) + geom_histogram(binwidth = 0.05) + labs(title="Distribution of sigmaPT",x="sigmaPT (median)") + coord_cartesian(ylim=c(0,175000),xlim=c(0,3))

ggplot(sigmaT_Cfull,aes(x=sigmaT_Cfull$p50)) + geom_histogram(binwidth = 0.05) + labs(title="Distribution of sigmaT",x="sigmaT (median)") + coord_cartesian(ylim=c(0,175000),xlim=c(0,3))

ggplot(sigmaE_Cfull,aes(x=sigmaE_Cfull$p50)) + geom_histogram(binwidth = 0.05) + labs(title="Distribution of sigmaE",x="sigmaE (median)") + coord_cartesian(ylim=c(0,175000),xlim=c(0,3))


### Log Ratios

ggplot(sigma,aes(x=log(sigma$sigmaP_Cfull.p50/sigma$sigmaPT_Cfull.p50))) + geom_histogram(binwidth = 0.2) + labs(title="Log Ratio of sigmaP/sigmaPT",x="log(sigmaP/sigmaPT)") + geom_vline(aes(xintercept=mean(log(sigma$sigmaP_Cfull.p50/sigma$sigmaPT_Cfull.p50))),col="blue",linetype="dashed") + annotate(geom="text",-4,70000,label=paste("Mean = ", round(mean(log(sigma$sigmaP_Cfull.p50/sigma$sigmaPT_Cfull.p50)),2)))

ggplot(sigma,aes(x=log(sigma$sigmaPT_Cfull.p50/sigma$sigmaT_Cfull.p50))) + geom_histogram(binwidth = 0.2) + labs(title="Log Ratio of sigmaPT/sigmaT",x="log(sigmaPT/sigmaT)") + geom_vline(aes(xintercept=mean(log(sigma$sigmaPT_Cfull.p50/sigma$sigmaT_Cfull.p50))),col="blue",linetype="dashed") + annotate(geom="text",-3,90000,label=paste("Mean = ", round(mean(log(sigma$sigmaPT_Cfull.p50/sigma$sigmaT_Cfull.p50)),2)))

ggplot(sigma,aes(x=log(sigma$sigmaP_Cfull.p50/sigma$sigmaT_Cfull.p50))) + geom_histogram(binwidth = 0.2) + labs(title="Log Ratio of sigmaP/sigmaT",x="log(sigmaP/sigmaT)") + geom_vline(aes(xintercept=mean(log(sigma$sigmaP_Cfull.p50/sigma$sigmaT_Cfull.p50))),col="blue",linetype="dashed") + annotate(geom="text",-4,70000,label=paste("Mean = ", round(mean(log(sigma$sigmaP_Cfull.p50/sigma$sigmaT_Cfull.p50)),2)))


### Correlation plots

ggplot(sigma,aes(x=sigma[,1],y=sigma[,2])) + geom_point(alpha = 0.01) + labs(x="sigmaP", y="sigmaPT", title="Correlation of sigmaP and sigmaPT") + theme(aspect.ratio = 1)
cor(sigma$sigmaP_Cfull.p50,sigma$sigmaPT_Cfull.p50)

ggplot(sigma,aes(x=sigma[,1],y=sigma[,3])) + geom_point(alpha = 0.01) + labs(x="sigmaP", y="sigmaT", title="Correlation of sigmaP and sigmaT") + coord_fixed(xlim=c(0,3),ylim=c(0,3)) + theme(aspect.ratio = 1)
cor(sigma$sigmaP_Cfull.p50,sigma$sigmaT_Cfull.p50)

ggplot(sigma,aes(x=sigma[,2],y=sigma[,3])) + geom_point(alpha = 0.01) + labs(x="sigmaPT", y="sigmaT", title="Correlation of sigmaPT and sigmaT") + coord_fixed(xlim=c(0,3),ylim=c(0,3)) + theme(aspect.ratio = 1)
cor(sigma$sigmaPT_Cfull.p50,sigma$sigmaT_Cfull.p50)



### Functional Region

load("myFA.Rdata")

PTRatio <- log(sigma$sigmaP_Cfull.p50/sigma$sigmaT_Cfull.p50)
PTTRatio <- log(sigma$sigmaPT_Cfull.p50/sigma$sigmaT_Cfull.p50)
PPTRatio <- log(sigma$sigmaP_Cfull.p50/sigma$sigmaPT_Cfull.p50)


# identify enhancers
Enhancer <- rep(NA,dim(FullAnnotation)[1])
for (i in 1:dim(FullAnnotation)[1]) {
  if (FullAnnotation$Phantom4_Enhancers[i] != '' | FullAnnotation$Phantom5_Enhancers[i] != '') {
    Enhancer[i] <- 1
  } else {
    Enhancer[i] <- 0
  }
}

# identify all possibilities in gene groups
temp <- strsplit(FullAnnotation$UCSC_RefGene_Group,";")
unique(unlist(temp))

# [1] "TSS1500" "Body"    "3'UTR"   "1stExon" "TSS200"  "5'UTR"   "5URT"    "3UTR"   
# [9] "ExonBnd"

Promoter1500 <- rep(NA,dim(FullAnnotation)[1])
Promoter200 <- rep(NA,dim(FullAnnotation)[1])
Body <- rep(NA,dim(FullAnnotation)[1])
Exon <- rep(NA,dim(FullAnnotation)[1])
UTR_5 <- rep(NA,dim(FullAnnotation)[1])
UTR_3 <- rep(NA,dim(FullAnnotation)[1])

for (i in 1:dim(FullAnnotation)[1]) {
  
 
  split <- strsplit(FullAnnotation$UCSC_RefGene_Group[i],";")
  
  if (sum(split[[1]] %in% "TSS1500") != 0) {
    
    Promoter1500[i] <- 1
    
  } else{
    Promoter1500[i] <- 0
  }
  
  if (sum(split[[1]] %in% "TSS200") != 0) {
    
    Promoter200[i] <- 1
    
  } else {
   
    Promoter200[i] <- 0
  }

  if (sum(split[[1]] %in% "Body") != 0 ) {
    Body[i] <- 1
  } else {
    Body[i] <- 0
  }

  if (sum(split[[1]] %in% "1stExon" | split[[1]] %in% "ExonBnd") != 0 ) {
    Exon[i] <- 1
  } else {
    Exon[i] <- 0
  }

  if (sum(split[[1]] %in% "5URT" | split[[1]] %in% "5'UTR") != 0 ) {
    UTR_5[i] <- 1
  } else {
    UTR_5[i] <- 0
  }

  if (sum(split[[1]] %in% "3UTR" | split[[1]] %in% "3'UTR") != 0 ) {
    UTR_3[i] <- 1
  } else {
    UTR_3[i] <- 0
  }

}

fun <- data.frame(Enhancer,Promoter1500,Promoter200,Body,Exon,UTR_3,UTR_5)

for (i in 1:dim(FullAnnotation)[1]) {
  if (sum(fun[i,]) == 0) {
    fun$Blank[i] <- 1
  } else {
    fun$Blank[i] <- 0
  }
}







hist(PTRatio,main="Distribution of Enhancers",
     xlab="log ratio of sigmaP/sigmaT",breaks=200,
     col=grey.colors(1,alpha=0.1))

hist(PTRatio[Enhancer==1],
     breaks=200,col="cyan", add=TRUE)

legend("topright",c("All Sites","Enhancer"),fill=c(grey.colors(1,alpha=0.1),"cyan"))
text(3,15000,labels=paste("Mean = ",round(mean(PTRatio[Enhancer==1]),2)))


hist(PTRatio,main="Distribution of TS1500",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=200,
     col=grey.colors(1,alpha=0.1))

hist(PTRatio[Promoter1500==1],
     breaks=200,col="coral", add=TRUE)

legend("topright",c("All Sites","TS1500"),fill=c(grey.colors(1,alpha=0.1),"coral"))
text(3,15000,labels=paste("Mean = ",round(mean(PTRatio[Promoter1500==1]),2)))


hist(PTRatio,main="Distribution of TS200",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=200,
     col=grey.colors(1,alpha=0.1))

hist(PTRatio[Promoter200==1],
     breaks=200,col="red", add=TRUE)

legend("topright",c("All Sites","TS200"),fill=c(grey.colors(1,alpha=0.1),"red"))
text(3,15000,labels=paste("Mean = ",round(mean(PTRatio[Promoter200==1]),2)))



hist(PTRatio,main="Distribution of Gene Body",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=200,
     col=grey.colors(1,alpha=0.1))

hist(PTRatio[Body==1],
     breaks=200,col="deeppink", add=TRUE)

legend("topright",c("All Sites","Gene Body"),fill=c(grey.colors(1,alpha=0.1),"deeppink"))
text(3,15000,labels=paste("Mean = ",round(mean(PTRatio[Body==1]),2)))


hist(PTRatio,main="Distribution of Exon Regions",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=200,
     col=grey.colors(1,alpha=0.1))

hist(PTRatio[Exon==1],
     breaks=200,col="darkorchid1", add=TRUE)

legend("topright",c("All Sites","Exon Region"),fill=c(grey.colors(1,alpha=0.1),"darkorchid1"))
text(3,15000,labels=paste("Mean = ",round(mean(PTRatio[Exon==1]),2)))




hist(PTRatio,main="Distribution of 3'UTR",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=200,
     col=grey.colors(1,alpha=0.1))

hist(PTRatio[UTR_3==1],
     breaks=200,col="dodgerblue", add=TRUE)

legend("topright",c("All Sites","3'UTR"),fill=c(grey.colors(1,alpha=0.1),"dodgerblue"))
text(3,15000,labels=paste("Mean = ",round(mean(PTRatio[UTR_3==1]),2)))


hist(PTRatio,main="Distribution of 5'UTR",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=200,
     col=grey.colors(1,alpha=0.1))

hist(PTRatio[UTR_5==1],
     breaks=200,col="gold", add=TRUE)

legend("topright",c("All Sites","5'UTR"),fill=c(grey.colors(1,alpha=0.1),"gold"))
text(3,15000,labels=paste("Mean = ",round(mean(PTRatio[UTR_5==1]),2)))



hist(PTRatio,main="Distribution of Blanks",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=200,
     col=grey.colors(1,alpha=0.1))

hist(PTRatio[fun$Blank==1],
     breaks=200,col="indianred1", add=TRUE)

legend("topright",c("All Sites","Blank"),fill=c(grey.colors(1,alpha=0.1),"indianred1"))
text(3,15000,labels=paste("Mean = ",round(mean(PTRatio[fun$Blank==1]),2)))



# Presentation used (only three regions)
hist(PTRatio,main="Distribution of Functional Regions",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=200,
     col=grey.colors(1,alpha=0.01))

hist(PTRatio[Promoter200==1],
     breaks=200,col=alpha("green",0.8), add=TRUE)

hist(PTRatio[Body==1],
     breaks=200,col=alpha("deeppink",0.5), add=TRUE)

hist(PTRatio[fun$Blank==1],
     breaks=200,col=alpha("yellow",0.3), add=TRUE)

legend("topright",c("All Sites","TSS200(Promoter)","Gene Body","Blank"),
       fill=c(grey.colors(1,alpha=0.1),alpha("green",0.8),alpha("deeppink",0.5),
              alpha("yellow",0.3)))

text(2.5,12500,labels=paste("Overall Mean = ",round(mean(PTRatio),2)))
text(2.5,11500,labels=paste("TSS200 Mean = ",round(mean(PTRatio[Promoter200==1]),2)))
text(2.5,10500,labels=paste("Body Mean = ",round(mean(PTRatio[Body==1]),2)))
text(2.5,9500,labels=paste("Blank Mean = ",round(mean(PTRatio[fun$Blank==1]),2)))


# Multiple Testing on Diff of Log Ratios for Functional Region Subgroups Within Gene

genelist <- c("APC","TP53","TTN","B2M","HLA-A","HLA-B")

for (i in genelist) {
  

tempPT <- PTRatio[geneInd(i)]
tempPPT <- PPTRatio[geneInd(i)]
tempPTT <- PTTRatio[geneInd(i)]

r1_enhancer <- tempPT[Enhancer[geneInd(i)]==1]
r2_enhancer <- tempPPT[Enhancer[geneInd(i)]==1]
r3_enhancer <- tempPTT[Enhancer[geneInd(i)]==1]

r1_200 <- tempPT[Promoter200[geneInd(i)]==1]
r2_200 <- tempPPT[Promoter200[geneInd(i)]==1]
r3_200 <- tempPTT[Promoter200[geneInd(i)]==1]

r1_1500 <- tempPT[Promoter1500[geneInd(i)]==1]
r2_1500 <- tempPPT[Promoter1500[geneInd(i)]==1]
r3_1500 <- tempPTT[Promoter1500[geneInd(i)]==1]

r1_exon <- tempPT[Exon[geneInd(i)]==1]
r2_exon <- tempPPT[Exon[geneInd(i)]==1]
r3_exon <- tempPTT[Exon[geneInd(i)]==1]

r1_body <- tempPT[Body[geneInd(i)]==1]
r2_body <- tempPPT[Body[geneInd(i)]==1]
r3_body <- tempPTT[Body[geneInd(i)]==1]

r1_3 <- tempPT[UTR_3[geneInd(i)]==1]
r2_3 <- tempPPT[UTR_3[geneInd(i)]==1]
r3_3 <- tempPTT[UTR_3[geneInd(i)]==1]

r1_5 <- tempPT[UTR_5[geneInd(i)]==1]
r2_5 <- tempPPT[UTR_5[geneInd(i)]==1]
r3_5 <- tempPTT[UTR_5[geneInd(i)]==1]

r1_b <- tempPT[fun$Blank[geneInd(i)]==1]
r2_b <- tempPPT[fun$Blank[geneInd(i)]==1]
r3_b <- tempPTT[fun$Blank[geneInd(i)]==1]


test.data1 <- data.frame(y=c(r1_enhancer,r1_1500,r1_200,r1_exon,r1_body,r1_3,r1_5,r1_b),
                        region=factor(rep(c("Enhancer","TS1500","TS200","Exon","Body",
                                            "3'UTR","5'UTR","Blank"),times=c(length(r1_enhancer),
                                                                             length(r1_1500),
                                                                     length(r1_200),
                                                                     length(r1_exon),
                                                                     length(r1_body),
                                                                     length(r1_3),
                                                                     length(r1_5),
                                                                     length(r1_b)))))


test.data2 <- data.frame(y=c(r2_enhancer,r2_1500,r2_200,r2_exon,r2_body,r2_3,r2_5,r2_b),
                         region=factor(rep(c("Enhancer","TS1500","TS200","Exon","Body",
                                             "3'UTR","5'UTR","Blank"),times=c(length(r2_enhancer),
                                                                              length(r2_1500),
                                                                              length(r2_200),
                                                                              length(r2_exon),
                                                                              length(r2_body),
                                                                              length(r2_3),
                                                                              length(r2_5),
                                                                              length(r2_b)))))

test.data3 <- data.frame(y=c(r3_enhancer,r3_1500,r3_200,r3_exon,r3_body,r3_3,r3_5,r3_b),
                         region=factor(rep(c("Enhancer","TS1500","TS200","Exon","Body",
                                             "3'UTR","5'UTR","Blank"),times=c(length(r3_enhancer),
                                                                              length(r3_1500),
                                                                              length(r3_200),
                                                                              length(r3_exon),
                                                                              length(r3_body),
                                                                              length(r3_3),
                                                                              length(r3_5),
                                                                              length(r3_b)))))

print(paste("Results of ",i))


# Dunn test for multiple comparison 
print("log P/T")
dunn.test(test.data1[,1],  test.data1[,2], method="bh")

print("log P/PT")
dunn.test(test.data2[,1],  test.data2[,2], method="bh")

print("log PT/T")
dunn.test(test.data3[,1],  test.data3[,2], method="bh")

}