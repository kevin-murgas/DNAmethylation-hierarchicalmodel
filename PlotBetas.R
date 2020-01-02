### plot beta values between specific samples

### INITIALIZE

library(ggplot2)
library(minfi)
library(MASS)

# here set the working directory that points to the data folder
# e.g. the folder with annotated data saved as "myFA.Rdata"
# please just comment out the directories already here, so each
# user can uncomment the setwd for them as we move code around

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# Yanlin's working directory
setwd("D:/DataPlus2017/Data")


##### LOAD DATA

# load fully annotated data (saved from LoadDataAndQC.R)
load("myFA.Rdata")


### mds plot

# mds plot
cutData <- FullAnnotation[-(1:8)]

pat <- data.frame(Sample.ID = colnames(cutData))

pat$Sample.ID = colnames(cutData)

aInd = c(1,2,14,15,21,22,29,30,52:55)
cInd = c(3:13,16:20,23:28,31,38,39,49,50,56:71)
oInd = c(32:37,40:48)

pat$type[aInd] <- "Adenomas"
pat$type[cInd] <- "Carcinomas"
pat$type[oInd] <- "Other"

patLabel <- pat$Sample.ID

cutData2 <- cutData[pat$type=="Adenomas"|pat$type=="Carcinomas"]

mdsPlot(as.matrix(cutData),numPositions = 1000000,sampNames=patLabel,sampGroups = pat$type,
        main = "PCA Plot Using Noob (with JA)", pch=16,cex=0.5)




##### COMPARE BETA

### tumor vs normal for patients H E K* W C J (have both)
# average tumor side A and B for tumor data

# H
betaHNor <- FullAnnotation$HN
betaHTum <- (FullAnnotation$HA + FullAnnotation$HB)/2
Hdata <- data.frame(x=betaHNor, y=betaHTum)
ggplot(Hdata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta H Normal", y="Beta H Tumor")
cor(betaHNor, betaHTum)
c(mean(betaHNor), mean(betaHTum), mean(betaHNor-betaHTum))

# E
betaENor <- FullAnnotation$EN
betaETum <- (FullAnnotation$EA + FullAnnotation$EB)/2
Edata <- data.frame(x=betaENor, y=betaETum)
ggplot(Edata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta E Normal", y="Beta E Tumor")
cor(betaENor, betaETum)
c(mean(betaENor), mean(betaETum), mean(betaENor-betaETum))

# K*
betaKnNor <- FullAnnotation$`KN (new)`
betaKnTum <- (FullAnnotation$`KA (new)` + FullAnnotation$`KB (new)`)/2
Kndata <- data.frame(x=betaKnNor, y=betaKnTum)
ggplot(Kndata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta K* Normal", y="Beta K* Tumor")
cor(betaKnNor, betaKnTum)
c(mean(betaKnNor), mean(betaKnTum), mean(betaKnNor-betaKnTum))

# W
betaWNor <- FullAnnotation$WN
betaWTum <- (FullAnnotation$WA + FullAnnotation$WB)/2
Wdata <- data.frame(x=betaWNor, y=betaWTum)
ggplot(Wdata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta W Normal", y="Beta W Tumor")
cor(betaWNor, betaWTum)
c(mean(betaWNor), mean(betaWTum), mean(betaWNor-betaWTum))

# C
betaCNor <- FullAnnotation$CN
betaCTum <- (FullAnnotation$CA + FullAnnotation$CB)/2
Cdata <- data.frame(x=betaCNor, y=betaCTum)
ggplot(Cdata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta C Normal", y="Beta C Tumor")
cor(betaCNor, betaCTum)
c(mean(betaCNor), mean(betaCTum), mean(betaCNor-betaCTum))

# J
betaJNor <- FullAnnotation$JN
betaJTum <- (FullAnnotation[,39] + FullAnnotation$JB)/2
Jdata <- data.frame(x=betaJNor, y=betaJTum)
ggplot(Jdata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta J Normal", y="Beta J Tumor")
cor(betaJNor, betaJTum)
c(mean(betaJNor), mean(betaJTum), mean(betaJNor-betaJTum))


### tumor side A vs side B

# K
betaKA <- FullAnnotation$KA
betaKB <- FullAnnotation$KB
Ktumor <- data.frame(x=betaKA, y=betaKB)
ggplot(Ktumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta K Side A", y="Beta K Side B")
cor(betaKA, betaKB)
c(mean(betaKA), mean(betaKB), mean(betaKA-betaKB))

# H
betaHA <- FullAnnotation$HA
betaHB <- FullAnnotation$HB
Htumor <- data.frame(x=betaHA, y=betaHB)
ggplot(Htumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta H Side A", y="Beta H Side B")
cor(betaHA, betaHB)

# E
betaEA <- FullAnnotation$EA
betaEB <- FullAnnotation$EB
Etumor <- data.frame(x=betaEA, y=betaEB)
ggplot(Etumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta E Side A", y="Beta E Side B")
cor(betaEA, betaEB)

# K*
betaKnA <- FullAnnotation$`KA (new)`
betaKnB <- FullAnnotation$`KB (new)`
Kntumor <- data.frame(x=betaKnA, y=betaKnB)
ggplot(Kntumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta K* Side A", y="Beta K* Side B")
cor(betaKnA, betaKnB)

# X
betaXA <- FullAnnotation$XA
betaXB <- FullAnnotation$XB
Xtumor <- data.frame(x=betaXA, y=betaXB)
ggplot(Xtumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta X Side A", y="Beta X Side B")
cor(betaXA, betaXB)

# W
betaWA <- FullAnnotation$WA
betaWB <- FullAnnotation$WB
Wtumor <- data.frame(x=betaWA, y=betaWB)
ggplot(Wtumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta W Side A", y="Beta W Side B")
cor(betaWA, betaWB)

# T
betaTA <- FullAnnotation$TA
betaTB <- FullAnnotation$TB
Ttumor <- data.frame(x=betaTA, y=betaTB)
ggplot(Ttumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta T Side A", y="Beta T Side B")
cor(betaTA, betaTB)

# S
betaSA <- FullAnnotation$SA
betaSB <- FullAnnotation$SB
Stumor <- data.frame(x=betaSA, y=betaSB)
ggplot(Stumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta S Side A", y="Beta S Side B")
cor(betaSA, betaSB)
c(mean(betaSA), mean(betaSB), mean(betaSA-betaSB))

# C
betaCA <- FullAnnotation$CA
betaCB <- FullAnnotation$CB
Ctumor <- data.frame(x=betaCA, y=betaCB)
ggplot(Ctumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta C Side A", y="Beta C Side B")
cor(betaCA, betaCB)
c(mean(betaCA), mean(betaCB), mean(betaCA-betaCB))

# M
betaMA <- FullAnnotation$MA
betaMB <- FullAnnotation$MB
Mtumor <- data.frame(x=betaMA, y=betaMB)
ggplot(Mtumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta M Side A", y="Beta M Side B")
cor(betaMA, betaMB)
c(mean(betaMA), mean(betaMB), mean(betaMA-betaMB))

# P
betaPA <- FullAnnotation$PA
betaPB <- FullAnnotation$PB
Ptumor <- data.frame(x=betaPA, y=betaPB)
ggplot(Ptumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta P Side A", y="Beta P Side B")
cor(betaPA, betaPB)
c(mean(betaPA), mean(betaPB), mean(betaPA-betaPB))

# J
betaJA <- FullAnnotation[,39]
betaJB <- FullAnnotation$JB
Jtumor <- data.frame(x=betaJA, y=betaJB)
ggplot(Jtumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta J Side A", y="Beta J Side B")
cor(betaJA, betaJB)
c(mean(betaJA), mean(betaJB), mean(betaJA-betaJB))

# O
betaOA <- FullAnnotation$OA
betaOB <- FullAnnotation$OB
Otumor <- data.frame(x=betaOA, y=betaOB)
ggplot(Otumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta O Side A", y="Beta O Side B")
cor(betaOA, betaOB)
c(mean(betaOA), mean(betaOB), mean(betaOA-betaOB))

# F
betaFA <- FullAnnotation$FA
betaFB <- FullAnnotation$FB
Ftumor <- data.frame(x=betaFA, y=betaFB)
ggplot(Ftumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta F Side A", y="Beta F Side B")
cor(betaFA, betaFB)
c(mean(betaFA), mean(betaFB), mean(betaFA-betaFB))

# D
betaDA <- FullAnnotation$DA
betaDB <- FullAnnotation$DB
Dtumor <- data.frame(x=betaDA, y=betaDB)
ggplot(Dtumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta D Side A", y="Beta D Side B")
cor(betaDA, betaDB)
c(mean(betaDA), mean(betaDB), mean(betaDA-betaDB))

# U
betaUA <- FullAnnotation$UA
betaUB <- FullAnnotation$UB
Utumor <- data.frame(x=betaUA, y=betaUB)
ggplot(Utumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta U Side A", y="Beta U Side B")
cor(betaUA, betaUB)
c(mean(betaUA), mean(betaUB), mean(betaUA-betaUB))


### Gland vs Gland

# KA
betaKAg1 <- FullAnnotation[,60]
betaKAg2 <- FullAnnotation[,61]
KAgland <- data.frame(x=betaKAg1, y=betaKAg2)
ggplot(KAgland,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta KA Gland 1", y="Beta KA Gland 2")
cor(betaKAg1, betaKAg2)
c(mean(betaKAg1), mean(betaKAg2), mean(betaKAg1-betaKAg2))

# KB
betaKBg1 <- FullAnnotation[,62]
betaKBg2 <- FullAnnotation[,63]
KBgland <- data.frame(x=betaKBg1, y=betaKBg2)
ggplot(KBgland,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta KB Gland 1", y="Beta KB Gland 2")
cor(betaKBg1, betaKBg2)
c(mean(betaKBg1), mean(betaKBg2), mean(betaKBg1-betaKBg2))

# CA
betaCAg1 <- FullAnnotation[,64]
betaCAg2 <- FullAnnotation[,65]
CAgland <- data.frame(x=betaCAg1, y=betaCAg2)
ggplot(CAgland,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta CA Gland 1", y="Beta CA Gland 2")
cor(betaCAg1, betaCAg2)
c(mean(betaCAg1), mean(betaCAg2), mean(betaCAg1-betaCAg2))

# CB
betaCBg1 <- FullAnnotation[,66]
betaCBg2 <- FullAnnotation[,67]
CBgland <- data.frame(x=betaCBg1, y=betaCBg2)
ggplot(CBgland,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta CB Gland 1", y="Beta CB Gland 2")
cor(betaCBg1, betaCBg2)
c(mean(betaCBg1), mean(betaCBg2), mean(betaCBg1-betaCBg2))

# FA
betaFAg1 <- FullAnnotation[,68]
betaFAg2 <- FullAnnotation[,69]
FAgland <- data.frame(x=betaFAg1, y=betaFAg2)
ggplot(FAgland,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta FA Gland 1", y="Beta FA Gland 2")
cor(betaFAg1, betaFAg2)
c(mean(betaFAg1), mean(betaFAg2), mean(betaFAg1-betaFAg2))

# FB
betaFBg1 <- FullAnnotation[,70]
betaFBg2 <- FullAnnotation[,71]
FBgland <- data.frame(x=betaFBg1, y=betaFBg2)
ggplot(FBgland,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta FB Gland 1", y="Beta FB Gland 2")
cor(betaFBg1, betaFBg2)
c(mean(betaFBg1), mean(betaFBg2), mean(betaFBg1-betaFBg2))

# HA
betaHAg1 <- FullAnnotation[,76]
betaHAg2 <- FullAnnotation[,77]
HAgland <- data.frame(x=betaHAg1, y=betaHAg2)
ggplot(HAgland,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta HA Gland 1", y="Beta HA Gland 2")
cor(betaHAg1, betaHAg2)
c(mean(betaHAg1), mean(betaHAg2), mean(betaHAg1-betaHAg2))

# HB
betaHBg1 <- FullAnnotation[,78]
betaHBg2 <- FullAnnotation[,79]
HBgland <- data.frame(x=betaHBg1, y=betaHBg2)
ggplot(HBgland,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta HB Gland 1", y="Beta HB Gland 2")
cor(betaHBg1, betaHBg2)
c(mean(betaHBg1), mean(betaHBg2), mean(betaHBg1-betaHBg2))


### Comparisons using only pat C and H
# C
betaCNor <- FullAnnotation$CN
betaCA <- FullAnnotation$CA
betaCB <- FullAnnotation$CB
betaCAg1 <- FullAnnotation[,64]
betaCAg2 <- FullAnnotation[,65]
betaCBg1 <- FullAnnotation[,66]
betaCBg2 <- FullAnnotation[,67]

# H
betaHNor <- FullAnnotation$HN
betaHA <- FullAnnotation$HA
betaHB <- FullAnnotation$HB
betaHAg1 <- FullAnnotation[,76]
betaHAg2 <- FullAnnotation[,77]
betaHBg1 <- FullAnnotation[,78]
betaHBg2 <- FullAnnotation[,79]


CvHnorm <- data.frame(x=betaCNor, y=betaHNor)
ggplot(CvHnorm,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta C Normal", y="Beta H Normal")
cor(betaCNor, betaHNor)

CvHA <- data.frame(x=betaCA, y=betaHA)
ggplot(CvHA,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta CA", y="Beta HA")
cor(betaCA, betaHA)

CvHB <- data.frame(x=betaCB, y=betaHB)
ggplot(CvHB,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta CB", y="Beta HB")
cor(betaCB, betaHB)

CvHAg1 <- data.frame(x=betaCAg1, y=betaHAg1)
ggplot(CvHAg1,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta CA gland 1", y="Beta HA gland 1")
cor(betaCAg1, betaHAg1)

CvHAg2 <- data.frame(x=betaCAg2, y=betaHAg2)
ggplot(CvHAg2,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta CA gland 2", y="Beta HA gland 2")
cor(betaCAg2, betaHAg2)

CvHBg1 <- data.frame(x=betaCBg1, y=betaHBg1)
ggplot(CvHBg1,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta CB gland 1", y="Beta HB gland 1")
cor(betaCBg1, betaHBg1)

CvHBg2 <- data.frame(x=betaCBg2, y=betaHBg2)
ggplot(CvHBg2,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta CB gland 2", y="Beta HB gland 2")
cor(betaCBg2, betaHBg2)


### CONTOUR PLOT
# right now this is not working for contour plot

z <- kde2d(HBgland[,1], HBgland[,2], h=0.5, n=30)
filled.contour(z, nlevels=15)

# C
betaCNor <- FullAnnotation$CN
betaCA <- FullAnnotation$CA
betaCB <- FullAnnotation$CB
betaCTum <- (betaCA+betaCB)/2
Ctissue <- data.frame(x=betaCNor, y=betaCTum)
zCtis <- kde2d(Ctissue[,1], Ctissue[,2], h=0.5, n=40)
filled.contour(zCtis, nlevels=20, xlab="Normal", ylab="Tumor")

Ctumor <- data.frame(x=betaCA, y=betaCB)
zCtum <- kde2d(Ctumor[,1], Ctumor[,2], h=0.5, n=40)
filled.contour(zCtum, zlim= ,nlevels=20, xlab="Bulk A", ylab="Bulk B")

cor(betaCNor, betaCTum)
c(mean(betaCNor), mean(betaCTum), mean(betaCNor-betaCTum))



ggplot(Ctumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta C Side A", y="Beta C Side B")
cor(betaCA, betaCB)
c(mean(betaCA), mean(betaCB), mean(betaCA-betaCB))


# H
betaHNor <- FullAnnotation$HN
betaHTum <- (FullAnnotation$HA + FullAnnotation$HB)/2
Hdata <- data.frame(x=betaHNor, y=betaHTum)
ggplot(Hdata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta H Normal", y="Beta H Tumor")
cor(betaHNor, betaHTum)
c(mean(betaHNor), mean(betaHTum), mean(betaHNor-betaHTum))

betaHA <- FullAnnotation$HA
betaHB <- FullAnnotation$HB
Htumor <- data.frame(x=betaHA, y=betaHB)
ggplot(Htumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta H Side A", y="Beta H Side B")
cor(betaHA, betaHB)