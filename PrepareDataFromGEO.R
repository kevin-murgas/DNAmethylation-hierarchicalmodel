# Script to load and prepare data for Hierarchical Methylation Conservation Model
# KA Murgas et al. 2021
# Data can be acquired from NCBI GEO dataset GSE166212
# data link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166212
# download the file "GSE166212_Normalized_Beta_values.xlsx"

# setwd to where the data file is located
setwd("/Users/kevinmurgas/Documents/Data+ project/GSE166212")
filename <- "GSE166212_Normalized_Beta_values.xlsx"

# use readxl library within tidyverse package (install if needed)
library(tidyverse)

# read first row to select relevant columns (skip detection pvals, only take beta_vals as numeric)
header <- read_xlsx(filename, n_max=1)
cols_skip <- if_else(str_detect(colnames(header), "Detection"), "skip", "numeric")
cols_skip[1] = "text" # first column is CpG site ID read as text

# read in data, skipping unnecessary columns, might take a minute
SampleData <- read_xlsx(filename, col_types = cols_skip)
samp_names <- colnames(SampleData) %>% str_split("_") %>% sapply("[[", 2)
samp_names[1] <- "IlmnID"
colnames(SampleData) <- samp_names
#save(SampleData, file="mySD.Rdata")

# read EPIC manifest annotation
# need to have downlaoded file "MethylationEPIC_v-1-0_B4.csv" (or more recent version) in the working directory
# available at: https://support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html
# these specified columns are to extract only relevant annotation information for each CpG site, more information is available in the manifest file
selectedCols <- c("character", rep("NULL", 10), "character", "integer", "NULL", rep("character", 2), "NULL", rep("character", 3), rep("NULL", 2), "character", rep("NULL", 24))
EPICchar <- read.csv("MethylationEPIC_v-1-0_B4.csv", skip=7) %>%
  select(IlmnID, CHR, MAPINFO, Strand, UCSC_RefGene_Name, UCSC_RefGene_Group, UCSC_CpG_Islands_Name, Relation_to_UCSC_CpG_Island, DMR)

# inner join EPICchar to SampleData
FullAnnotation = inner_join(EPICchar, SampleData, by="IlmnID") %>% arrange(IlmnID)

# save output file as "myFA_bulkonly.Rdata"
save(FullAnnotation, file = "myFA_bulkonly.Rdata")
