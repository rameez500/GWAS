############################################
##### Pre-Imputation Processing and QC #####
############################################

## Note: Anytime you see data, data, or n, those should be specific to your data you're working on

## Analyst:
## Name of data:
## Location:
## Genome build:
## Platform:
## Strand:
## Number variants:
## Number of samples:
## Original filename:

## Load libraries needed in this pipeline
library(qqman)
library(data.table)
library(ggplot2)
library(readxl)
library(tidyr)
library(CMplot)
library(car)
library(corrplot)
library(ggcorrplot)
library(Hmisc)
library(car)
library(foreach)



## Take a look at the original .bim and .fam files
system("wc -l filepath-to-original-data/data.bim") # n-variants

## Look at markers in data -- check to see which build they correspond to using dbSNP: https://www.ncbi.nlm.nih.gov/snp/
bim <- fread("filepath-to-original-data/data.bim")
bim[3000:3010,] ## Better to look in the middle of the file instead of just the begining/head 
## Paste bim sample here

## Run the .bim file through Chipendium to determine the microarray platform: http://mccarthy.well.ox.ac.uk/chipendium/ui/

## Look at the structure of the fam file -- will use this information after imputation to correct IIDs
fam <- fread("filepath-to-original-data/data.fam")
head(fam)
## Paste fam sample here

nrow(fam) # n-people


###################
##### PART A ######
###################

#####  Step 0: Lift data to NCBI37 #####
#####  Step 0: Lift data to NCBI37 #####

## Download strand files from Raynor website: http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/index.html
## Strand file:

## Required parameters: update_build.sh <bed-file-stem> <strand-file> <output-file-stem>
##  1. The original bed stem (not including file extension suffix)
##  2. The strand file to apply
##  3. The new stem for output

system("sh /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/update_build.sh original_data PATH-TO-STRAND-FILE data_37")

system("wc -l data_37.fam") # n-people
system("wc -l data_37.bim") # n-variants

## Note: If you do not have to lift to NCBI37, run original file through plink and generate new .bed/.bim/.fam files with _37 appended to filename
## Change --bfile to --file if not in plink binary format 
system("plink --bfile filepath-to-original-data/data --make-bed --out data_37")


##### Step 1: QC the Sample data #####
##### Step 1: QC the Sample data #####

## Screen for genotyping rate >0.05
system("plink --bfile data_37 --geno 0.05 --threads 20 --make-bed --out data_37_geno")
#n variants and n people pass filters and QC

system("wc -l data_37_geno.fam") # n-people
system("wc -l data_37_geno.bim") # n-variants

## Screen for MAF > 0.1
system("plink --bfile data_37_geno --maf 0.1 --threads 20 --make-bed --out data_37_maf")
#n variants and n people pass filters and QC

system("wc -l data_37_maf.fam") # n-people
system("wc -l data_37_maf.bim") # n-variants


##### Step 2: Merge Reference and Sample data to identify overlapping SNPs #####
##### Step 2: Merge Reference and Sample data to identify overlapping SNPs #####

## 1000 Genomes LD pruned data (50 10 .5) file (original file from SVS)
library(data.table)
kg_bim <- fread("/projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/IMPUTE2_1KG_Reference_All_Biallelic_SNPS.bim")
dim(kg_bim) # 2050463       13
head(kg_bim)
kg_bim_snps <- as.matrix(kg_bim$V2)

## Read in your sample data .bim
data_bim <- fread("data_37_maf.bim")
head(data_bim)
nrow(data_bim) # n variants
data_bim_snps <- as.matrix(data_bim$V2)

## Determine overlap between 1KG and sample data
table(kg_bim_snps %in% data_bim_snps)
#    FALSE    TRUE
# n  n
overlap <- kg_bim_snps[kg_bim_snps %in% data_bim_snps]
length(overlap) # n
library(data.table)
fwrite(as.list(overlap), "data_updated_1kg_ref_snps.txt", sep = " ", quote = FALSE, col.names = FALSE)


##### Step 3: Restrict datas to overlapping SNPs #####
##### Step 3: Restrict datas to overlapping SNPs #####

## Create new binary 1KG reference files of only overlapping SNPs
system("plink --bfile /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/IMPUTE2_1KG_Reference_All_Biallelic_SNPS --extract data_updated_1kg_ref_snps.txt  --make-bed --out 1kg_reference_for_data")
#n variants and n people pass filters and QC
system("wc -l 1kg_reference_for_data.bim") # n

## Trim down the data to match reference sample
system("plink --bfile data_37_maf --extract data_updated_1kg_ref_snps.txt  --make-bed --out data_match_1kgRef")
#n variants and n people pass filters and QC

system("wc -l data_match_1kgRef.fam") # n-people
system("wc -l data_match_1kgRef.bim") # n-variants


##### Step 4: Rayner Tool #####
##### Step 4: Rayner Tool #####

## Make frequency file
system("plink --bfile data_match_1kgRef --freq --out data_match_1kgRef")

## Run tool
system("perl /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC-1000G-check-bim-v4.2.11_Oct2019/HRC-1000G-check-bim.pl -b data_match_1kgRef.bim -f data_match_1kgRef.frq -r /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/1000GP_Phase3_combined.legend -g -p ALL")
# Options Set:
# Reference Panel:             1000G
# Bim filename:                data_match_1kgRef.bim
# Reference filename:          /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/1000GP_Phase3_combined.legend
# Allele frequencies filename: data_match_1kgRef.frq
# Population for 1000G:        ALL

system("mv Run-plink.sh Run_plink_data.sh")
system("chmod u+x Run_plink_data.sh")
system("./Run_plink_data.sh")

## Check how many are excluded
system("wc -l Exclude-data_match_1kgRef-1000G.txt") # n are excluded
system("wc -l data_match_1kgRef-updated.bim") # n are retained
system("wc -l data_match_1kgRef-updated.fam") # n are retained


##### Step 5: Prepare file for PCA #####
##### Step 5: Prepare file for PCA #####

## Create new binary 1KG reference files of only overlapping SNPs
system("plink --bfile 1kg_reference_for_data --exclude Exclude-data_match_1kgRef-1000G.txt  --make-bed --out data_match_1kgRef_raynor")
#n variants and n people pass filters and QC


## Merge updated data with 1KG for PCA
system("plink --bfile data_match_1kgRef-updated --bmerge data_match_1kgRef_raynor --threads 20 --make-bed --out data_match_1kgRef-updated_1KG")
#n variants and n people pass filters and QC

system("wc -l data_match_1kgRef-updated_1KG.bim") # n-variants
system("wc -l data_match_1kgRef-updated_1KG.fam") # n-people + 1kg


##### Step 6: FlashPCA #####
##### Step 6: FlashPCA #####

## Run PCA
system("flashpca --bfile data_match_1kgRef-updated_1KG --ndim 20 --suffix _data.txt --numthreads 20 > pca_data.log &")

## Get the PCs
pcs <-fread("pcs_data.txt")
dim(pcs) # n-people+1kg    22

## Read in 1KG Super Populations data for plotting
thousand_index <- fread("/projects/bga_lab/DATA_REPOSITORIES/1000Genome/1000Genome_labels.csv")
sub_thous <- thousand_index[,c(1,3)]
colnames(sub_thous)[1] <- c("IID")
table(sub_thous$GROUP, useNA = "ifany")
# AFR AMR EAS EUR SAS
# 661 347 504 503 489

## Merge Super Population labels into PC data
PC_20_merged <- merge(pcs, sub_thous, by = "IID", all.x= TRUE ) #all.x=TRUE allows you to keep the sample with NAs for Population Group
table(PC_20_merged$GROUP, useNA = "ifany") #should now have NAs from the c1 sample
# AFR  AMR  EAS  EUR  SAS <NA>
# 661  347  504  503  489 n-people 
PC_20_merged$GROUP[is.na(PC_20_merged$GROUP)] <- "DATA" #rename the NAs to  data
table(PC_20_merged$GROUP, useNA = "ifany")
# AFR    AMR    EAS    EUR DATA    SAS
# 661    347    504    503   n-people    489

## Plot Ancestry based on 1KG reference groups
library(ggplot2)
library(tidyverse)

ggplot(PC_20_merged, aes(x=PC1, y=PC2, color=GROUP)) +
  geom_point(size=0.2) +
  theme_classic()
ggsave("data_pc2vpc1.png", width = 5, height = 5, units = "in")

ggplot(PC_20_merged, aes(x=PC2, y=PC3, color=GROUP)) +
  geom_point(size=0.2) +
  theme_classic()
ggsave("data_pc2vpc3.png", width = 5, height = 5, units = "in")

ggplot(PC_20_merged, aes(x=PC3, y=PC4, color=GROUP)) +
  geom_point(size=0.2) +
  theme_classic()
ggsave("data_pc3vpc4.png", width = 5, height = 5, units = "in")


##### Step 7: Select Ancestral Groups ####
##### Step 7: Select Ancestral Groups ####

## Select out DATA participants
data <- PC_20_merged[(PC_20_merged$GROUP == "DATA"), ]
head(data)
dim(data) # n     23

## Calculate means, standard deviations for each Super Population
## AFR
just_afr <- PC_20_merged[(PC_20_merged$GROUP == "AFR"), ]
head(just_afr)
dim(just_afr) # n     23
afr_m1 <- mean(just_afr$PC1)
afr_m2 <- mean(just_afr$PC2)
afr_m3 <- mean(just_afr$PC3)
afr_sd1 <-sd(just_afr$PC1)
afr_sd2 <-sd(just_afr$PC2)
afr_sd3 <-sd(just_afr$PC3)

afr_means <- rbind(afr_m1, afr_m2, afr_m3)
afr_sds <- rbind(afr_sd1, afr_sd2, afr_sd3)
afr_stats <- cbind(afr_means, afr_sds)
colnames(afr_stats) <- c("mean", "stand_dev")
afr_stats
#       mean   stand_dev
# afr_m1  0.xx 0.xx
# afr_m2  0.xx 0.xx
# afr_m3 -0.xx 0.xx

## EUR
just_eur <- PC_20_merged[(PC_20_merged$GROUP == "EUR"), ]
head(just_eur)
dim(just_eur) # 503     23
eur_m1 <- mean(just_eur$PC1)
eur_m2 <- mean(just_eur$PC2)
eur_m3 <- mean(just_eur$PC3)
eur_sd1 <-sd(just_eur$PC1)
eur_sd2 <-sd(just_eur$PC2)
eur_sd3 <-sd(just_eur$PC3)

eur_means <- rbind(eur_m1, eur_m2, eur_m3)
eur_sds <- rbind(eur_sd1, eur_sd2, eur_sd3)
eur_stats <- cbind(eur_means, eur_sds)
colnames(eur_stats) <- c("mean", "stand_dev")
eur_stats
#              mean   stand_dev
# eur_m1 -0.xx 0.xx
# eur_m2  0.xx 0.xx
# eur_m3 -0.xx 0.xx


## Compute the thresholds for outliers +/- 2 SD's away from mean
## AFR
afr_UP_thresh_1 <- (afr_stats[1,1] + (2*afr_stats[1,2]))
afr_LO_thresh_1 <- (afr_stats[1,1] - (2*afr_stats[1,2]))
afr_UP_thresh_2 <- (afr_stats[2,1] + (2*afr_stats[2,2]))
afr_LO_thresh_2 <- (afr_stats[2,1] - (2*afr_stats[2,2]))
afr_UP_thresh_3 <- (afr_stats[3,1] + (2*afr_stats[3,2]))
afr_LO_thresh_3 <- (afr_stats[3,1] - (2*afr_stats[3,2]))

afr_thresholds <- rbind(afr_UP_thresh_1, afr_LO_thresh_1, afr_UP_thresh_2, afr_LO_thresh_2,afr_UP_thresh_3, afr_LO_thresh_3)
afr_thresholds
# afr_UP_thresh_1  0.xx
# afr_LO_thresh_1  0.xx
# afr_UP_thresh_2  0.xx
# afr_LO_thresh_2 -0.xx
# afr_UP_thresh_3  0.xx
# afr_LO_thresh_3 -0.xx

## EUR
eur_UP_thresh_1 <- (eur_stats[1,1] + (2*eur_stats[1,2]))
eur_LO_thresh_1 <- (eur_stats[1,1] - (2*eur_stats[1,2]))
eur_UP_thresh_2 <- (eur_stats[2,1] + (2*eur_stats[2,2]))
eur_LO_thresh_2 <- (eur_stats[2,1] - (2*eur_stats[2,2]))
eur_UP_thresh_3 <- (eur_stats[3,1] + (2*eur_stats[3,2]))
eur_LO_thresh_3 <- (eur_stats[3,1] - (2*eur_stats[3,2]))

eur_thresholds <- rbind(eur_UP_thresh_1, eur_LO_thresh_1, eur_UP_thresh_2, eur_LO_thresh_2,eur_UP_thresh_3, eur_LO_thresh_3)
eur_thresholds
# eur_UP_thresh_1 -0.xx
# eur_LO_thresh_1 -0.xx
# eur_UP_thresh_2  0.xx
# eur_LO_thresh_2  0.xx
# eur_UP_thresh_3  0.xx
# eur_LO_thresh_3 -0.xx

all_thresholds <- cbind(afr_thresholds, eur_thresholds)
colnames(all_thresholds) <- c("afr", "eur")
rownames(all_thresholds) <- c("upper_thresh_pc1", "lower_thresh_pc1", "upper_thresh_pc2", "lower_thresh_pc2", "upper_thresh_pc3", "lower_thresh_pc3")
all_thresholds

## Remove Outliers for each Super Population of interest: AFR, EUR
## AFR
dim(data) # n     23
data_afr_1up <- data[data$PC1 <= all_thresholds[1,1],]
dim(data_afr_1up) # n   23
head(data_afr_1up)
data_afr_1lo <- data_afr_1up[data_afr_1up$PC1 >= all_thresholds[2,1],]
dim(data_afr_1lo) # n   23
data_afr_2up <- data_afr_1lo[data_afr_1lo$PC2 <= all_thresholds[3,1],]
dim(data_afr_2up) # n   23
data_afr_2lo <- data_afr_2up[data_afr_2up$PC2 >= all_thresholds[4,1],]
dim(data_afr_2lo) # n   23
data_afr_3up <- data_afr_2lo[data_afr_2lo$PC3 <= all_thresholds[5,1],]
dim(data_afr_3up) # n   23
data_afr_3lo <- data_afr_3up[data_afr_3up$PC3 >= all_thresholds[6,1],]
dim(data_afr_3lo) # n   23


## EUR
dim(data) # n     23
data_eur_1up <- data[data$PC1 <= all_thresholds[1,2],]
dim(data_eur_1up) # n     23
head(data_eur_1up)
data_eur_1lo <- data_eur_1up[data_eur_1up$PC1 >= all_thresholds[2,2],]
dim(data_eur_1lo) # n     23
data_eur_2up <- data_eur_1lo[data_eur_1lo$PC2 <= all_thresholds[3,2],]
dim(data_eur_2up) # n     23
data_eur_2lo <- data_eur_2up[data_eur_2up$PC2 >= all_thresholds[4,2],]
dim(data_eur_2lo) # n     23
data_eur_3up <- data_eur_2lo[data_eur_2lo$PC3 <= all_thresholds[5,2],]
dim(data_eur_3up) # n     23
data_eur_3lo <- data_eur_3up[data_eur_3up$PC3 >= all_thresholds[6,2],]
dim(data_eur_3lo) # n     23

## Write out subject keep lists for each group of interest
data_afr_2lo_iids <- data_afr_2lo[,c(2,1)] # n
nrow(data_afr_2lo_iids)
data_afr_3lo_iids <- data_afr_3lo[,c(2,1)] # n
nrow(data_afr_3lo_iids)
data_eur_2lo_iids <- data_eur_2lo[,c(2,1)] # n
nrow(data_eur_2lo_iids)
data_eur_3lo_iids <- data_eur_3lo[,c(2,1)] # n
nrow(data_eur_3lo_iids)

## Save them out
fwrite(data_afr_2lo_iids, "data_afr_2pcs_keeplist.txt", col.names = FALSE, sep = " ")
fwrite(data_afr_3lo_iids, "data_afr_3pcs_keeplist.txt", col.names = FALSE, sep = " ")
fwrite(data_eur_2lo_iids, "data_eur_2pcs_keeplist.txt", col.names = FALSE, sep = " ")
fwrite(data_eur_3lo_iids, "data_eur_3pcs_keeplist.txt", col.names = FALSE, sep = " ")

## Get PCA files restricted to AFR subjects
afr_pcs <- data[data$IID %in% data_afr_3lo_iids$IID,]
dim(afr_pcs) # n     23

## Get PCA files restricted to EUR subjects
eur_pcs <- data[data$IID %in% data_eur_3lo_iids$IID,]
dim(eur_pcs) # n     23

## Write out
fwrite(afr_pcs, "data_afr_pcs.txt")
fwrite(eur_pcs, "data_eur_pcs.txt")

## Implement multidimensional scaling as done in SVS
# see the website below for more information on calculation used
# http://doc.goldenhelix.com/SVS/latest/svsmanual/numeric_data_quality.html#multoutlier


## Calculate medians and quartiles for each column
# Afr
#  Calculate quartile values for each PC (convert to matrix so will rbind correctly with other stats)
afr_quart_pc1 <- as.matrix(quantile(afr_pcs$PC1))
afr_quart_pc2 <- as.matrix(quantile(afr_pcs$PC2))
afr_quart_pc3 <- as.matrix(quantile(afr_pcs$PC3))

#interquartile ranges for each PC
afr_iqr_pc1 <- afr_quart_pc1[4,1] - afr_quart_pc1[2,1] #Q3 - Q1
afr_iqr_pc2 <- afr_quart_pc2[4,1] - afr_quart_pc2[2,1] #Q3 - Q1
afr_iqr_pc3 <- afr_quart_pc3[4,1] - afr_quart_pc3[2,1] #Q3 - Q1

afr_pc1_stats <- rbind(afr_iqr_pc1, afr_quart_pc1)
afr_pc2_stats <- rbind(afr_iqr_pc2, afr_quart_pc2)
afr_pc3_stats <- rbind(afr_iqr_pc3, afr_quart_pc3)

afr_stats2 <- cbind(afr_pc1_stats, afr_pc2_stats, afr_pc3_stats)
afr_stats2
rownames(afr_stats2) <- c("iqr", "min", "q1", "median", "q3", "max")
colnames(afr_stats2) <- c("PC1", "PC2", "PC3")
afr_stats2

#                PC1         PC2          PC3
# iqr    
# min    
# q1     
# median 
# q3     
# max  

## Calculate medians and quartiles for each column
# EUR
#  Calculate quartile values for each PC (convert to matrix so will rbind correctly with other stats)
eur_quart_pc1 <- as.matrix(quantile(eur_pcs$PC1))
eur_quart_pc2 <- as.matrix(quantile(eur_pcs$PC2))
eur_quart_pc3 <- as.matrix(quantile(eur_pcs$PC3))

#interquartile ranges for each PC
eur_iqr_pc1 <- eur_quart_pc1[4,1] - eur_quart_pc1[2,1] #Q3 - Q1
eur_iqr_pc2 <- eur_quart_pc2[4,1] - eur_quart_pc2[2,1] #Q3 - Q1
eur_iqr_pc3 <- eur_quart_pc3[4,1] - eur_quart_pc3[2,1] #Q3 - Q1

eur_pc1_stats <- rbind(eur_iqr_pc1, eur_quart_pc1)
eur_pc2_stats <- rbind(eur_iqr_pc2, eur_quart_pc2)
eur_pc3_stats <- rbind(eur_iqr_pc3, eur_quart_pc3)

eur_stats2 <- cbind(eur_pc1_stats, eur_pc2_stats, eur_pc3_stats)
eur_stats2
rownames(eur_stats2) <- c("iqr", "min", "q1", "median", "q3", "max")
colnames(eur_stats2) <- c("PC1", "PC2", "PC3")
eur_stats2

#                PC1         PC2          PC3
# iqr    
# min    
# q1     
# median 
# q3     
# max 

###  Multidimensional Outlier Detection  ###
## AFR
## Step A - Calculate the distance for each sample (row) of a pc.txt
## square root (sum (value - median)^2)
dim(afr_pcs) # n   23
head(afr_pcs)

m_afr_pcs <- afr_pcs[,1:5]
head(m__pcs)
m_afr_pcs$distance <- 999 # This is just a place holder value, will be replaced with actual distance in for loop below
m_afr_pcs$outlier <- 999 # This is another place holder
for (i in 1:nrow(m_afr_pcs)){
  m_afr_pcs[i,6] <- sqrt(((m_afr_pcs[i,3]-afr_stats2[4,1])^2)+((m_afr_pcs[i,4]-afr_stats2[4,2])^2)+((m_afr_pcs[i,5]-afr_stats2[4,3])^2))
}
head(m_afr_pcs)


## Step B - calculate threshold; Outlier >= threshold
threshold_afr <- (sqrt((afr_stats2[5,1]^2)+(afr_stats2[5,2]^2)+(afr_stats2[5,3]^2)) + 1.5*(sqrt((afr_stats2[1,1]^2) + (afr_stats2[1,2]^2) + (afr_stats2[1,3]^2))))
threshold_afr # 0.xx

## Step C - classify as outlier
for (i in 1:nrow(m_afr_pcs)){
  m_afr_pcs[i,7] <- ifelse(m_afr_pcs[i,6]>= threshold_afr,1,0)
}
head(m_afr_pcs)

table(m_afr_pcs$outlier) # 0 outliers
#    0
# n

## Double check that we don't have outliers
max(m_afr_pcs$distance) # 0.xx (this is less than the threshold)

## Select out individuals to keep
head(m_afr_pcs)
afr_mds_keep_data <- m_afr_pcs[m_afr_pcs$outlier == '0',]
nrow(afr_mds_keep_data)
fwrite(afr_mds_keep_data[,c(2,1)], "data_afr_3pcs_keeplist_mds.txt", col.names = FALSE, sep = " ", row.names = F)


###  Multidimensional Outlier Detection  ###
## EUR
## Step A - Calculate the distance for each sample (row) of a pc.txt
## square root (sum (value - median)^2)
dim(eur_pcs) # n   23
head(eur_pcs)

m_eur_pcs <- eur_pcs[,1:5]
head(m_eur_pcs)
m_eur_pcs$distance <- 999 # This is just a place holder value, will be replaced with actual distance in for loop below
m_eur_pcs$outlier <- 999 # This is another place holder
for (i in 1:nrow(m_eur_pcs)){
  m_eur_pcs[i,6] <- sqrt(((m_eur_pcs[i,3]-eur_stats2[4,1])^2)+((m_eur_pcs[i,4]-eur_stats2[4,2])^2)+((m_eur_pcs[i,5]-eur_stats2[4,3])^2))
}
head(m_eur_pcs)


## Step B - calculate threshold; Outlier >= threshold
threshold_eur <- (sqrt((eur_stats2[5,1]^2)+(eur_stats2[5,2]^2)+(eur_stats2[5,3]^2)) + 1.5*(sqrt((eur_stats2[1,1]^2) + (eur_stats2[1,2]^2) + (eur_stats2[1,3]^2))))
threshold_eur # 0.xx

## Step C - classify as outlier
for (i in 1:nrow(m_eur_pcs)){
  m_eur_pcs[i,7] <- ifelse(m_eur_pcs[i,6]>= threshold_eur,1,0)
}
head(m_eur_pcs)

table(m_eur_pcs$outlier) # 0 outliers
#    0
# n

## Double check that we don't have outliers
max(m_eur_pcs$distance) # 0.xx (this is less than the threshold)

## Select out individuals to keep
head(m_eur_pcs)
eur_mds_keep_data <- m_eur_pcs[m_eur_pcs$outlier == '0',]
nrow(eur_mds_keep_data)
fwrite(eur_mds_keep_data[,c(2,1)], "data_eur_3pcs_keeplist_mds.txt", col.names = FALSE, sep = " ", row.names = F)


#####  Step 8: Extract super populations from sample data #####
#####  Step 8: Extract super populations from sample data #####
## AFR
## Based on the individuals identified in Part A, create subsets for each super population (the file created after Step 0)
system("plink --bfile data_37 --keep data_afr_3pcs_keeplist_mds.txt --make-bed --out data_37_afr")
# n variants and n people pass filters and QC

system("wc -l data_37_afr.fam")
system("wc -l data_37_afr.bim")

## EUR
## Based on the individuals identified in Part A, create subsets for each super population (the file created after Step 0)
system("plink --bfile data_37 --keep data_eur_3pcs_keeplist_mds.txt --make-bed --out data_37_eur")
# n variants and n people pass filters and QC

system("wc -l data_37_eur.fam")
system("wc -l data_37_eur.bim")

##### Step 9: Rayner Tool #####
##### Step 9: Rayner Tool #####

## AFR
## Make frequency file
system("plink --bfile data_37_afr --freq --out data_37_afr")

## Run pipeline
system("perl /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC-1000G-check-bim-v4.2.11_Oct2019/HRC-1000G-check-bim.pl -b data_37_afr.bim -f data_37_afr.frq -r /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/CAAPA_FILES/all.caapa.sorted.txt -h -t .1")

# Options Set:
# Reference Panel:             HRC (note: the software prints out HRC, but we are using the CAAPA file)
# Bim filename:                data_37_afr.bim
# Reference filename:          /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/CAAPA_FILES/all.caapa.sorted.txt
# Allele frequencies filename: data_37_afr.frq
# Chromosome flag set:         No
# Allele frequency threshold:  0.1

system("mv Run-plink.sh Run_plink_data_37_afr.sh")
system("chmod u+x Run_plink_data_37_afr.sh")
system("./Run_plink_data_37_afr.sh")

## Check how many are excluded
system("wc -l Exclude-data_37_afr-HRC.txt") # n are excluded
system("wc -l data_37_afr-updated.bim") # n are retained
system("wc -l data_37_afr-updated.fam") # n are retained

## EUR
## Make frequency file
system("plink --bfile data_37_eur --freq --out data_37_eur")

## Run pipeline
system("perl /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC-1000G-check-bim-v4.2.11_Oct2019/HRC-1000G-check-bim.pl -b data_37_eur.bim -f data_37_eur.frq -r /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC_FILES/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -t 0.1")

# Options Set:
# Reference Panel:             HRC 
# Bim filename:                data_37_eur.bim
# Reference filename:          /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC.r1-1.GRCh37.wgs.mac5.sites.tab
# Allele frequencies filename: data_37_eur.frq
# Chromosome flag set:         No
# Allele frequency threshold:  0.1

system("mv Run-plink.sh Run_plink_data_37_eur.sh")
system("chmod u+x Run_plink_data_37_eur.sh")
system("./Run_plink_data_37_eur.sh")

## Check how many are excluded
system("wc -l Exclude-data_37_eur-HRC.txt") # n are excluded
system("wc -l data_37_eur-updated.bim") # n are retained
system("wc -l data_37_eur-updated.fam") # n are retained



##### Step 10: QC the Sample data by chromosome #####
##### Step 10: QC the Sample data by chromosome #####

## Create list of chrm files to use
AFR_chrms <- foreach(i=1:22) %do%paste("data_37_afr-updated-chr", i, sep = "")
EUR_chrms <- foreach(i=1:22) %do%paste("data_37_eur-updated-chr", i, sep = "")

## Loop through each chrm and QC individually
## AFR ##
## Screen for genotyping rate >0.05
foreach(i=AFR_chrms) %do% {
  cmd_geno <- paste0("plink --bfile ", i, " --geno 0.05 --threads 20 --make-bed --out ", i, "_geno", sep = "")
  system(cmd_geno)}

# Get total number of markers across all files after --geno
markers_geno <- foreach(i=AFR_chrms) %do% {nrow(fread(paste(i, "_geno.bim", sep = "")))} # n's by chrm   
Reduce("+", markers_geno) ## Total n
# n markers

# Get total number of people across all files after --geno
unique(unlist(foreach(i=AFR_chrms) %do% {nrow(fread(paste(i, "_geno.fam", sep = "")))}))
# n people

## Screen for MAF > 0.01
foreach(i=AFR_chrms) %do% {
  cmd_maf <- paste0("plink --bfile ", i, "_geno --maf 0.01 --threads 20 --make-bed --out ", i, "_maf", sep = "")
  system(cmd_maf)}

# Get total number of markers across all files after --maf
markers_maf <- foreach(i=AFR_chrms) %do% {nrow(fread(paste(i, "_maf.bim", sep = "")))}    
Reduce("+",markers_maf) ## Total n
# n markers

# Get total number of people across all files after --maf
unique(unlist(foreach(i=AFR_chrms) %do% {nrow(fread(paste(i, "_maf.fam", sep = "")))}))
# n people


## Screen for missingness by person
foreach(i=AFR_chrms) %do% {
  cmd_maf <- paste0("plink --bfile ", i, "_maf --mind 0.1 --threads 20 --make-bed --out ", i, "_mind", sep = "")
  system(cmd_maf)}

# Get total number of markers across all files after --maf
markers_mind <- foreach(i=AFR_chrms) %do% {nrow(fread(paste(i, "_mind.bim", sep = "")))}    
Reduce("+",markers_mind)## Total n
# n markers

# Get total number of people across all files after --maf
unique(unlist(foreach(i=AFR_chrms) %do% {nrow(fread(paste(i, "_mind.fam", sep = "")))}))
# n people


## EUR
## Screen for genotyping rate >0.05
foreach(i=EUR_chrms) %do% {
  cmd_geno <- paste0("plink --bfile ", i, " --geno 0.05 --threads 20 --make-bed --out ", i, "_geno", sep = "")
  system(cmd_geno)}

# Get total number of markers across all files after --geno
markers_geno <- foreach(i=EUR_chrms) %do% {nrow(fread(paste(i, "_geno.bim", sep = "")))}    
Reduce("+", markers_geno)## Total n
# n markers

# Get total number of people across all files after --geno
unique(unlist(foreach(i=EUR_chrms) %do% {nrow(fread(paste(i, "_geno.fam", sep = "")))}))
# n people

## Screen for MAF > 0.01
foreach(i=EUR_chrms) %do% {
  cmd_maf <- paste0("plink --bfile ", i, "_geno --maf 0.01 --threads 20 --make-bed --out ", i, "_maf", sep = "")
  system(cmd_maf)}

# Get total number of markers across all files after --maf
markers_maf <- foreach(i=EUR_chrms) %do% {nrow(fread(paste(i, "_maf.bim", sep = "")))}    
Reduce("+",markers_maf)## Total n
# n markers

# Get total number of people across all files after --maf
unique(unlist(foreach(i=EUR_chrms) %do% {nrow(fread(paste(i, "_maf.fam", sep = "")))}))
# n people

## Screen for missingness by person
foreach(i=EUR_chrms) %do% {
  cmd_maf <- paste0("plink --bfile ", i, "_maf --mind 0.1 --threads 20 --make-bed --out ", i, "_mind", sep = "")
  system(cmd_maf)}

# Get total number of markers across all files after --mind
markers_mind <- foreach(i=EUR_chrms) %do% {nrow(fread(paste(i, "_mind.bim", sep = "")))}    
Reduce("+",markers_mind)## Total n
# n markers

# Get total number of people across all files after --mind
unique(unlist(foreach(i=EUR_chrms) %do% {nrow(fread(paste(i, "_mind.fam", sep = "")))}))
# n people


##### Step 11: Prep data for imputation #####
##### Step 11: Prep data for imputation #####

## AFR
## Convert to vcf (may take a few moments)
chrm <- 1:22
x2 <- NULL
for(i in 1:22) {
  x2 <- rbind(x2,paste("plink --bfile data_37_afr-updated-chr",chrm[i],"_mind --recode vcf --out data_37_afr-updated-chr",chrm[i],sep=""))
}
write.table(x2,file="data_37_afr-updated_bychrm_tovcf.sh",row.names=FALSE,col.names=FALSE,quote=FALSE)
system("chmod +x data_37_afr-updated_bychrm_tovcf.sh")
system("nohup ./data_37_afr-updated_bychrm_tovcf.sh > data_37_afr-updated_bychrm_tovcf.log &")

##### TOPMED IMPUTATION SERVER #####
####If the vcf file has to be lifted to version 38, use the following command to update the CHROM column in VCF file before submitting to the TOPMed imputation sever
#Run the following command in terminal as a bash command
#for i in {1..22}; do  echo '$i chr$i' >> chr_name_conv.txt; bcftools annotate --rename-chrs chr_name_conv.txt data_37_ancestry-updated-chr$i.vcf -o data_37_ancestry-updated-withCHR-chr$i.vcf; done
#Replace "data_37_afr-updated-chr" to "data_37_afr-updated-withCHR-chr" in the next step if you are recoding CHROM column

## Create sorted vcf (may take a few moments)
x3 <- NULL
for(i in 1:22) {
  x3 <- rbind(x3,paste("vcf-sort data_37_afr-updated-chr",chrm[i],".vcf | bgzip -c > data_37_afr-updated-chr",chrm[i],".vcf.gz",sep=""))
}
write.table(x3,file="data_37_afr-updated_tosortedvcf.sh",row.names=FALSE,col.names=FALSE,quote=FALSE)
system("chmod +x data_37_afr-updated_tosortedvcf.sh")
system("nohup ./data_37_afr-updated_tosortedvcf.sh > data_37_afr-updated_tosortedvcf.log &")

## EUR
## Convert to vcf (may take a few moments)
chrm <- 1:22
x2 <- NULL
for(i in 1:22) {
  x2 <- rbind(x2,paste("plink --bfile data_37_eur-updated-chr",chrm[i],"_mind --recode vcf --out data_37_eur-updated-chr",chrm[i],sep=""))
}
write.table(x2,file="data_37_eur-updated_bychrm_tovcf.sh",row.names=FALSE,col.names=FALSE,quote=FALSE)
system("chmod +x data_37_eur-updated_bychrm_tovcf.sh")
system("nohup ./data_37_eur-updated_bychrm_tovcf.sh > data_37_eur-updated_bychrm_tovcf.log &")

##### TOPMED IMPUTATION SERVER #####
####If the vcf file has to be lifted to version 38, use the following command to update the CHROM column in VCF file before submitting to the TOPMed imputation sever
#Run the following command in terminal as a bash command
#for i in {1..22}; do  echo '$i chr$i' >> chr_name_conv.txt; bcftools annotate --rename-chrs chr_name_conv.txt data_37_ancestry-updated-chr$i.vcf -o data_37_ancestry-updated-withCHR-chr$i.vcf; done
#Replace "data_37_eur-updated-chr" to "data_37_eur-updated-withCHR-chr" in the next step if you are recoding CHROM column

## Create sorted vcf (may take a few moments)
x3 <- NULL
for(i in 1:22) {
  x3 <- rbind(x3,paste("vcf-sort data_37_eur-updated-chr",chrm[i],".vcf | bgzip -c > data_37_eur-updated-chr",chrm[i],".vcf.gz",sep=""))
}
write.table(x3,file="data_37_eur-updated_tosortedvcf.sh",row.names=FALSE,col.names=FALSE,quote=FALSE)
system("chmod +x data_37_eur-updated_tosortedvcf.sh")
system("nohup ./data_37_eur-updated_tosortedvcf.sh > data_37_eur-updated_tosortedvcf.log &")


###########################
##### READY TO IMPUTE #####
###########################

# Server version:

### Copy and Paste Information from Imputation Server:

## AFR ##
## AFR ##

# 22 valid VCF file(s) found.
#
# Samples: n
# Chromosomes: 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9
# SNPs: n
# Chunks: n
# Datatype: unphased
# Reference Panel: hrc.r1.1.2016
# Phasing: shapeit


## EUR ##
## EUR ##

# 22 valid VCF file(s) found.
#
# Samples: n
# Chromosomes: 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9
# SNPs: n
# Chunks: n
# Datatype: unphased
# Reference Panel: hrc.r1.1.2016
# Phasing: shapeit