#################################################################################
######## This program is used for reading in dataset and qualify control ########
######### Refer to the RMarkdown file (coming soon) for final results ###########
#################################################################################


library(minfi)
library(minfiData)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(gdata) #gdata package needs Perl language installed
library(RColorBrewer)

# here set the working directory that points to the data folder
# e.g. the folder with all datasets in it, should contain all the
# 1337-1387 folders
# please just comment out the directories already here, so each
# user can uncomment the setwd for them as we move code around

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# Yanlin's working directory
setwd("D:/DataPlus2017/Data")

############### Function used for read in data ####################

# base_dir = directory containing basenames file
# targets = name of the dataset after reading in basenames file
# pat_dir = directory containing patient information file
# pat_file = WHOLE name (including directory) of the dataset after reading in patient information file
# work_name = name of the RGChannelSet

read.fun <- function(base_dir,pat_file) {
  # read in patient data from the .xlsx file
  pat_name <- read.xls(pat_file)
  
  # Extract targets
  targets_name <- data.frame(pat_name[,"Complete.Barcode"])
  colnames(targets_name) <- "Basement"
  
  # Read in all .idat files
  targets_name$Basement <- file.path(base_dir, targets_name$Basement)
  work_name <- read.metharray(targets_name$Basement, verbose = TRUE)
  
  # Push the RGChannelSet to the global environment
  return(work_name)
}


############################### Read in data  #################################

## 1337 ##
base_dir_1337 <- "1337_Shibata EPIC DNA methylation data package/IDAT FILES"
pat_file_1337 <- "1337_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1337 (Shibata-8).xls"
work_1337 <- read.fun(base_dir_1337,pat_file_1337)

## 1345 ##
base_dir_1345 <- "1345_Shibata EPIC DNA methylation data package/IDAT FILES"
pat_file_1345 <- "1345_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1345 (Shibata-16).xlsx"
work_1345 <- read.fun(base_dir_1345,pat_file_1345)

# Two beadchips in plate 1345, read in seperately
pat_1345 <- read.xls(pat_file_1345)
targets_1345 <- data.frame(pat_1345[,"Complete.Barcode"])
colnames(targets_1345) <- 'Basement'

targets_1345$Basement <- file.path(base_dir_1345,targets_1345$Basement)
# Beadchip 1
work_1345a <- read.metharray(targets_1345$Basement[1:8],verbose=TRUE)
# Beadchip 2
work_1345b <- read.metharray(targets_1345$Basement[9:16],verbose = TRUE)

## 1350 ##
base_dir_1350 <- "1350_SHIBATA EPIC DNA METHYLATION DATA PACKAGE/IDAT files"
pat_file_1350 <- "1350_SHIBATA EPIC DNA METHYLATION DATA PACKAGE/SAMPLE-ARRAY MAPPING/1350 (Shibata-8).xlsx"
work_1350 <- read.fun(base_dir_1350,pat_file_1350)
 
## 1357 ##
base_dir_1357 <- "1357_Shibata EPIC DNA methylation data package/IDAT FILES"
pat_file_1357 <- "1357_Shibata EPIC DNA Methylation Data Package/SAMPLE-ARRAY MAPPING/1357 (Shibata-8).xlsx"
work_1357 <- read.fun(base_dir_1357,pat_file_1357)

## 1360 ##
base_dir_1360 <- "1360_Shibata EPIC Data Package/IDAT FILES"
pat_file_1360 <- "1360_Shibata EPIC Data Package/SAMPLE-ARRAY MAPPING/1360 (Shibata-8).xlsx"
work_1360 <- read.fun(base_dir_1360,pat_file_1360)

## 1378 ##
base_dir_1378 <- "1378_Shibata EPIC Data Package/IDAT FILES"
pat_file_1378 <- "1378_Shibata EPIC Data Package/SAMPLE-ARRAY MAPPING/1378 (Shbata-8).xls"
work_1378 <- read.fun(base_dir_1378,pat_file_1378)

## 1385 ##
base_dir_1385 <- "1385_Shibata EPIC Data Package/IDAT FILES"
pat_file_1385 <- "1385_Shibata EPIC Data Package/SAMPLE-ARRAY MAPPING/1385 (Shibata-8).xlsx"
work_1385 <- read.fun(base_dir_1385,pat_file_1385)
 
## 1387 ##
base_dir_1387 <- "1387_Shibata EPIC DNA methylation data package/IDAT FILES"
pat_file_1387 <- "1387_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1387 (Shibata-8).xls"
work_1387 <- read.fun(base_dir_1387,pat_file_1387)

## 1586 ##
base_dir_1586 <- "1586_Shibata EPIC DNA methylation data package/IDATs"
pat_file_1586 <- "1586_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1586 (Shibata-16 FFPE).xlsx"
#work_1586 <- read.fun(base_dir_1586,pat_file_1586)

# Two beadchips in plate 1586, read in seperately
pat_1586 <- read.xls(pat_file_1586)
targets_1586 <- data.frame(pat_1586[,"Complete.Barcode"])
colnames(targets_1586) <- 'Basement'

targets_1586$Basement <- file.path(base_dir_1586,targets_1586$Basement)
# Beadchip 1
work_1586 <- read.metharray(targets_1586$Basement[1:8],verbose=TRUE)
# Beadchip 2 is glands, leave out for now

## 1639 ##
base_dir_1639 <- "1639_Shibata EPIC DNA methylation data package/IDAT FILES"
pat_file_1639 <- "1639_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1639 (Shibata-8).xlsx"
work_1639 <- read.fun(base_dir_1639,pat_file_1639)

################################## Bad Sample in 1337 ####################################

# Convert to MethylSet
mset_1337 <- preprocessRaw(work_1337)

# Check the methylation and unmethylation signals
head(getMeth(mset_1337))
head(getUnmeth((mset_1337)))

# Convert to RatioSet containing beta values and Mvalues(log(M/U))
rset_1337 <- ratioConvert(mset_1337,what = "both")

# Get beta
beta_1337 <- getBeta(rset_1337)

# Check beta values
head(beta_1337)

#200360140022_R07C01 (JB) has bad quality
dim(beta_1337)[1] == sum(is.na(beta_1337[,7]))

### read in new data without JB ###
base_dir_1337_new <- "1337_Shibata EPIC DNA methylation data package/IDAT FILES"
targets_1337_new <- read.csv("1337_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1337_Shibata-targets.csv", as.is = TRUE)
targets_1337_new$Basename <- file.path(base_dir_1337_new, targets_1337_new$Basename)
work_1337_new <- read.metharray(targets_1337_new$Basename, verbose = TRUE)


############################## Preprocessing/Normalization ###############################

### 1337 ###
## noob ##
noob_1337 <- preprocessNoob(work_1337_new)
rset_1337_noob <- ratioConvert(noob_1337,what="both")
beta_1337_noob <- getBeta(rset_1337_noob)
## Genome Studio ##
lumi_1337 <- preprocessIllumina(work_1337_new,bg.correct = TRUE,
                                normalize = "controls",reference = 2)
rset_1337_lumi <- ratioConvert(lumi_1337,what="both")
beta_1337_lumi <- getBeta(rset_1337_lumi)


### 1345 ###
## noob ##
noob_1345 <- preprocessNoob(work_1345)
noob_1345a <- preprocessNoob(work_1345a)
noob_1345b <- preprocessNoob(work_1345b)
## Genome Studio ##
lumi_1345 <- preprocessIllumina(work_1345,bg.correct = TRUE,
                                normalize = "controls",reference = 2)
lumi_1345a <- preprocessIllumina(work_1345a,bg.correct = TRUE,
                                normalize = "controls",reference = 2)
lumi_1345b <- preprocessIllumina(work_1345b,bg.correct = TRUE,
                                normalize = "controls",reference = 2)

### 1350 ###
## noob ##
noob_1350 <- preprocessNoob(work_1350)

## Genome Studio ##
lumi_1350 <- preprocessIllumina(work_1350,bg.correct = TRUE,
                                normalize = "controls",reference = 2)

### 1357 ###
## noob ##
noob_1357 <- preprocessNoob(work_1357)
## Genome Studio ##
lumi_1357 <- preprocessIllumina(work_1357,bg.correct = TRUE,
                                normalize = "controls",reference = 2)

### 1360 ###
## noob ##
noob_1360 <- preprocessNoob(work_1360)
## Genome Studio ##
lumi_1360 <- preprocessIllumina(work_1360,bg.correct = TRUE,
                                normalize = "controls",reference = 2)

### 1378 ###
## noob ##
noob_1378 <- preprocessNoob(work_1378)
## Genome Studio ##
lumi_1378 <- preprocessIllumina(work_1378,bg.correct = TRUE,
                                normalize = "controls",reference = 2)

### 1385 ###
## noob ##
noob_1385 <- preprocessNoob(work_1385)
## Genome Studio ##
lumi_1385 <- preprocessIllumina(work_1385,bg.correct = TRUE,
                                normalize = "controls",reference = 2)

### 1387 ###
## noob ##
noob_1387 <- preprocessNoob(work_1387)
## Genome Studio ##
lumi_1387 <- preprocessIllumina(work_1387,bg.correct = TRUE,
                                normalize = "controls",reference = 2)

### 1586 ###
## noob ##
noob_1586 <- preprocessNoob(work_1586a)
## Genome Studio ##
lumi_1586 <- preprocessIllumina(work_1586a,bg.correct = TRUE,
                                normalize = "controls",reference = 2)

### 1639 ###
## noob ##
noob_1639 <- preprocessNoob(work_1639)
## Genome Studio ##
lumi_1639 <- preprocessIllumina(work_1639,bg.correct = TRUE,
                                normalize = "controls",reference = 2)

#################################### Quality Control #####################################

####### Density Plots #########

## 1337 ##
pat_1337 <- read.xls(pat_file_1337)
pat_1337 <- pat_1337[-7,] # remove bad sample JB

par(mfrow=c(1,2),mar=c(4,2,4,1))
densityPlot(noob_1337,sampGroups = pat_1337$Sample_No,
            main = "Corrected Beta Values Using 'noob'",legend = FALSE)
legend("topright",legend = levels(pat_1337$Sample_No),
       text.col = brewer.pal(8,"Dark2"),cex=0.5)
densityPlot(lumi_1337,sampGroups = pat_1337$Sample_No,
            main = "Corrected Beta Values Using 'Illumina'",legend = FALSE)
legend("topright",legend = levels(pat_1337$Sample_No),
       text.col = brewer.pal(8,"Dark2"),cex=0.5)

## 1345 ##
pat_1345 <- read.xls(pat_file_1345)

par(mfrow=c(1,2),mar=c(4,2,4,1))
densityPlot(noob_1345,sampGroups = pat_1345$SAMPLE.ID,
            main = "Corrected Beta Values Using 'noob'",legend = FALSE)
legend("top",legend = levels(pat_1345$SAMPLE.ID),
       text.col = brewer.pal(8,"Dark2"),cex=0.5)
densityPlot(lumi_1345,sampGroups = pat_1345$SAMPLE.ID,
            main = "Corrected Beta Values Using 'Illumina'",legend = FALSE)
legend("top",legend = levels(pat_1345$SAMPLE.ID),
       text.col = brewer.pal(8,"Dark2"),cex=0.5)

## 1350 ##
pat_1350 <- read.xls(pat_file_1350)

par(mfrow=c(1,2),mar=c(4,2,4,1))
densityPlot(noob_1350,sampGroups = pat_1350$SAMPLE.ID,
            main = "Corrected Beta Values Using 'noob'",legend = FALSE)
legend("topright",legend = levels(pat_1350$SAMPLE.ID),
       text.col = brewer.pal(8,"Dark2"),cex=0.7)
densityPlot(lumi_1350,sampGroups = pat_1350$SAMPLE.ID,
            main = "Corrected Beta Values Using 'Illumina'",legend = FALSE)
legend("topright",legend = levels(pat_1350$SAMPLE.ID),
       text.col = brewer.pal(8,"Dark2"),cex=0.7)

## 1357 ## 
pat_1357 <- read.xls(pat_file_1357)

par(mfrow=c(1,2),mar=c(4,2,4,1))
densityPlot(noob_1357,sampGroups = pat_1357$Sample.ID,
            main = "Corrected Beta Values Using 'noob'",legend = FALSE)
legend("top",legend = levels(pat_1357$Sample.ID),
       text.col = brewer.pal(8,"Dark2"),cex=0.7)
densityPlot(lumi_1357,sampGroups = pat_1357$Sample.ID,
            main = "Corrected Beta Values Using 'Illumina'",legend = FALSE)
legend("top",legend = levels(pat_1357$Sample.ID),
       text.col = brewer.pal(8,"Dark2"),cex=0.7)

## 1360 ##
pat_1360 <- read.xls(pat_file_1360)

par(mfrow=c(1,2),mar=c(4,2,4,1))
densityPlot(noob_1360,sampGroups = pat_1360$Sample.ID,
            main = "Corrected Beta Values Using 'noob'",legend = FALSE)
legend("top",legend = levels(pat_1360$Sample.ID),
       text.col = brewer.pal(8,"Dark2"),cex=0.7)
densityPlot(lumi_1360,sampGroups = pat_1360$Sample.ID,
            main = "Corrected Beta Values Using 'Illumina'",legend = FALSE)
legend("top",legend = levels(pat_1360$Sample.ID),
       text.col = brewer.pal(8,"Dark2"),cex=0.7)

## 1378 ##
pat_1378 <- read.xls(pat_file_1378)

par(mfrow=c(1,2),mar=c(4,2,4,1))
densityPlot(noob_1378,sampGroups = pat_1378$Tube.Label,
            main = "Corrected Beta Values Using 'noob'",legend = FALSE)
legend("top",legend = levels(pat_1378$Tube.Label),
       text.col = brewer.pal(8,"Dark2"),cex=0.7)
densityPlot(lumi_1378,sampGroups = pat_1378$Tube.Label,
            main = "Corrected Beta Values Using 'Illumina'",legend = FALSE)
legend("top",legend = levels(pat_1378$Tube.Label),
       text.col = brewer.pal(8,"Dark2"),cex=0.7)

## 1385 ##
pat_1385 <- read.xls(pat_file_1385)

par(mfrow=c(1,2),mar=c(4,2,4,1))
densityPlot(noob_1385,sampGroups = pat_1385$Tube.Label,
            main = "Corrected Beta Values Using 'noob'",legend = FALSE)
legend("top",legend = levels(pat_1385$Tube.Label),
       text.col = brewer.pal(8,"Dark2"),cex=0.7)
densityPlot(lumi_1385,sampGroups = pat_1385$Tube.Label,
            main = "Corrected Beta Values Using 'Illumina'",legend = FALSE)
legend("top",legend = levels(pat_1385$Tube.Label),
       text.col = brewer.pal(8,"Dark2"),cex=0.7)

## 1387 ##
pat_1387 <- read.xls(pat_file_1387)

par(mfrow=c(1,2),mar=c(4,2,4,1))
densityPlot(noob_1387,sampGroups = pat_1387$Sample.ID,
            main = "Corrected Beta Values Using 'noob'",legend = FALSE)
legend("top",legend = levels(pat_1387$Sample.ID),
       text.col = brewer.pal(8,"Dark2"),cex=0.7)
densityPlot(lumi_1387,sampGroups = pat_1387$Sample.ID,
            main = "Corrected Beta Values Using 'Illumina'",legend = FALSE)
legend("top",legend = levels(pat_1387$Sample.ID),
       text.col = brewer.pal(8,"Dark2"),cex=0.7)

## 1586 ##
pat_1586 <- read.xls(pat_file_1586)
pat_1586 <- droplevels(pat_1586[1:8,]) # only use first 8

par(mfrow=c(1,2),mar=c(4,2,4,1))
densityPlot(noob_1586,sampGroups = pat_1586$USC.ID,
            main = "Corrected Beta Values Using 'noob'",legend = FALSE)
legend("topright",legend = levels(pat_1586$USC.ID),
       text.col = brewer.pal(8,"Dark2"),cex=0.5)
densityPlot(lumi_1586,sampGroups = pat_1586$USC.ID,
            main = "Corrected Beta Values Using 'Illumina'",legend = FALSE)
legend("topright",legend = levels(pat_1586$USC.ID),
       text.col = brewer.pal(8,"Dark2"),cex=0.5)

## 1639 ##
pat_1639 <- read.xls(pat_file_1639)

# because a patient in 1639 is named "NA", it loaded as <NA> so manually change to "NA"
levels(pat_1639$Tube.Label)<-c(levels(pat_1639$Tube.Label),"NA") # add "NA" level to factor first
pat_1639$Tube.Label[which(is.na(pat_1639$Tube.Label))]<-"NA"

par(mfrow=c(1,2),mar=c(4,2,4,1))
densityPlot(noob_1639,sampGroups = pat_1639$Tube.Label,
            main = "Corrected Beta Values Using 'noob'",legend = FALSE)
legend("topright",legend = levels(pat_1639$Tube.Label),
       text.col = brewer.pal(8,"Dark2"),cex=0.5)
densityPlot(lumi_1639,sampGroups = pat_1639$Tube.Label,
            main = "Corrected Beta Values Using 'Illumina'",legend = FALSE)
legend("topright",legend = levels(pat_1639$Tube.Label),
       text.col = brewer.pal(8,"Dark2"),cex=0.5)

######### Median Plot ############

qc_1337_noob <- getQC(noob_1337)
qc_1345_noob <- getQC(noob_1345)
qc_1350_noob <- getQC(noob_1350)
qc_1357_noob <- getQC(noob_1357)
qc_1360_noob <- getQC(noob_1360)
qc_1378_noob <- getQC(noob_1378)
qc_1385_noob <- getQC(noob_1385)
qc_1387_noob <- getQC(noob_1387)
qc_1586_noob <- getQC(noob_1586)
qc_1639_noob <- getQC(noob_1639)

qc_noob <- rbind(qc_1337_noob,qc_1345_noob,qc_1350_noob,qc_1357_noob,qc_1360_noob,
                 qc_1378_noob,qc_1385_noob,qc_1387_noob,qc_1586_noob,qc_1639_noob)

par(mfrow=c(1,1))
plotQC(qc_noob, main="dd")

qc_1337_lumi <- getQC(lumi_1337)
qc_1345_lumi <- getQC(lumi_1345)
qc_1350_lumi <- getQC(lumi_1350)
qc_1357_lumi <- getQC(lumi_1357)
qc_1360_lumi <- getQC(lumi_1360)
qc_1378_lumi <- getQC(lumi_1378)
qc_1385_lumi <- getQC(lumi_1385)
qc_1387_lumi <- getQC(lumi_1387)
qc_1586_lumi <- getQC(lumi_1586)
qc_1639_lumi <- getQC(lumi_1639)

qc_lumi <- rbind(qc_1337_lumi,qc_1345_lumi,qc_1350_lumi,qc_1357_lumi,qc_1360_lumi,
                 qc_1378_lumi,qc_1385_lumi,qc_1387_lumi,qc_1586_lumi,qc_1639_lumi)

par(mfrow=c(1,1))
plotQC(qc_lumi)


# Plot sex
gmset_1337 <- mapToGenome(noob_1337)
sex_1337 <- getSex(gmset_1337)
plotSex(sex_1337)

# MDS plot
mdsPlot(noob_1337)

all_noob <- NULL
all_noob <- combineArrays(noob_1337,noob_1345,noob_1350,noob_1357,noob_1360,
                          noob_1378,noob_1385,noob_1387,noob_1586,noob_1639)

# Percent of probes
per_fun <- function (work) {
  mset <- preprocessRaw(work)
  rset <- ratioConvert(mset,what = "both")
  beta <- getBeta(rset)
  per <- 1-apply(is.na(beta),2,sum)/dim(beta)[1]
  return(per)
}

per_1337 <- per_fun(work_1337) # work_1337_new
per_1345 <- per_fun(work_1345)
per_1350 <- per_fun(work_1350)
per_1357 <- per_fun(work_1357)
per_1360 <- per_fun(work_1360)
per_1378 <- per_fun(work_1378)
per_1385 <- per_fun(work_1385)
per_1387 <- per_fun(work_1387)
per_1586 <- per_fun(work_1586)
per_1639 <- per_fun(work_1639)

#################################### Annotation #####################################

# here we select, from the EPIC characterizing data, the columns for:
#                 ILmnID,                       CHR,         MAPINFO(position), Strand, UCSC_RefGene_Name,   USCS_RefGene_Group, CpG Island name, Relation_to_UCSC_CpG_Island, DMR
selectedCols <- c("character", rep("NULL", 10), "character", "integer", "NULL", rep("character", 2), "NULL", rep("character", 3), rep("NULL", 2), "character", rep("NULL", 24))
EPICchar <- read.csv("EPIC MANIFEST AND SUPPORTING INFORMATION/MethylationEPIC_v-1-0_B2.csv", as.is=TRUE, colClasses = selectedCols)

## then add each sample
# data from noob preprocessed data
SampleData <- data.frame(getBeta(noob_1337), getBeta(noob_1345), getBeta(noob_1350),
                         getBeta(noob_1357), getBeta(noob_1360), getBeta(noob_1378),
                         getBeta(noob_1385), getBeta(noob_1387), getBeta(noob_1586),
                         getBeta(noob_1639))
# rewrite sample names from patient file
sampleNames <- unlist(list(pat_1337$Sample_No, pat_1345$SAMPLE.ID, pat_1350$SAMPLE.ID,
                           pat_1357$Sample.ID, pat_1360$Sample.ID, pat_1378$Tube.Label,
                           pat_1385$Tube.Label, pat_1387$Sample.ID, pat_1586$USC.ID, pat_1639$Tube.Label))
colnames(SampleData) <- sampleNames[-7] # store sample names, removing bad sample JA in 1337
# reorder sample data by cpg site (row name)
SampleData <- SampleData[order(row.names(SampleData)), ]

# only choose rows with cpg sites in preprocessed data
DataChar <- EPICchar[EPICchar$IlmnID %in% row.names(SampleData),]

# save/load to save time
save(DataChar,file="myDC.Rdata")
save(SampleData,file="mySD.Rdata")
load("myDC.Rdata")
load("mySD.Rdata")

# merge two, then sort by IlmnID
FullAnnotation <- cbind(DataChar, SampleData)

save(FullAnnotation,file="myFA.Rdata")
load("myFA.Rdata")

########## Bulk-Samples Only (hard-coded) ##########
# extract only bulk tumors (manually, only works with original FA and 2 new data sets)
bulkInds <- c(1:14,16:40,47:48,58:59,73:76,87:96)
FullAnnotation <- FullAnnotation[,bulkInds]

save(FullAnnotation,file="myFA_bulkonly.Rdata")
load("myFA_bulkonly.Rdata")


######### Make a separate annotation using minfi's mapToGenome
library(minfi)
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

## 1639 ## Take 1639 and convert raw data to MethylSet
base_dir_1639 <- "1639_Shibata EPIC DNA methylation data package/IDAT FILES"
pat_file_1639 <- "1639_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1639 (Shibata-8).xlsx"
work_1639 <- read.fun(base_dir_1639,pat_file_1639)
raw_1639 <- preprocessRaw(work_1639)

map <- mapToGenome(raw_1639)
gr <- granges(map, use.names=TRUE)
anno <- getAnnotation(map)

# DMP finder (can compare to mu/betaT later)
temp <- as.matrix(FullAnnotation[,10:57])
rownames(temp) <- FullAnnotation$IlmnID
patientLabel <- substr(colnames(temp),1,1)
patientLabel[10:12] <- "K*"
sideLabel <- substr(colnames(temp),2,2)
tissueLabel <- sideLabel
tissueLabel[!(sideLabel %in% c("N"))] <- "T"   #replace A/B/D/M with T
tumorIndicator <- 1*(tissueLabel=="T")
dmp <- dmpFinder(temp, pheno=tissueLabel, type = "categorical")