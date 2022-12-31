#############################################
##### Post-Imputation Processing and QC #####  
#############################################

## analyst: 
## date:
## Name of dataset:
## Location of imputed data
## Ancestry you're working with: 
## change: "data" to your dataset and "ancestry" to AFR or EUR (or admixed)

## Load libraries needed in this script
library(data.table)
library(dplyr)
library(tidyr)

### Information from Imputation Server: job-xxxx ####

# Statistics:
#   Alternative allele frequency > 0.5 sites: 0
# Reference Overlap: 100.00 %
# Match: 392,499
# Allele switch: 165,887
# Strand flip: 0
# Strand flip and allele switch: 0
# A/T, C/G genotypes: 0
# Filtered sites:
#   Filter flag set: 0
# Invalid alleles: 0
# Multiallelic sites: 0
# Duplicated sites: 0
# NonSNP sites: 0
# Monomorphic sites: 0
# Allele mismatch: 664
# SNPs call rate < 90%: 0
# 
# Excluded sites in total: 664
# Remaining sites in total: 558,386
# See [NOT AVAILABLE] for details
# Typed only sites: 18
# See [NOT AVAILABLE] for details
# 
# Warning: 1 Chunk(s) excluded: < 3 SNPs (see [NOT AVAILABLE] for details).
# Remaining chunk(s): 153

###################
##### PART G ######
###################

## Data located: /projects/bga_lab/DATA_REPOSITORIES/dbGaP/ACTIVE/ADDHealth/IMPUTATION/POSTIMP/

##### Step 0: Download data from Michigan Server #####
##### Step 0: Download data from Michigan Server #####

## Imputation password: password

##### Step 1: Unzip each chr folder #####
##### Step 1: Unzip each chr folder #####

chrm <- 1:22
x <- NULL
for(i in 1:22) {
  x <- rbind(x,paste("7z x chr_",chrm[i],".zip -p'password'",sep=""))
  
}
write.table(x,file="unzip_data.sh",row.names=FALSE,col.names=FALSE,quote=FALSE)
system("dos2unix unzip_data.sh")
system("chmod +x unzip_data.sh")
system("nohup ./unzip_data.sh > unzip_data.log &")

##### Step 2: Unzip each .info file #####
##### Step 2: Unzip each .info file #####

x2 <- NULL
for(i in 1:22) {
  x2 <- rbind(x2,paste("gunzip chr",chrm[i],".info.gz",sep=""))
}
write.table(x2,file="unzip_info_data.sh",row.names=FALSE,col.names=FALSE,quote=FALSE)
system("dos2unix unzip_info_data.sh")
system("chmod +x unzip_info_data.sh")
system("nohup ./unzip_info_data.sh > unzip_info_data.log &")


##### Step 3: Prepare imputed files for merging #####
##### Step 3: Prepare imputed files for merging #####

## Read in .info file to create lists of markers to screen
## Merge the info files together

# Read in all chrm info files
chr1 <- fread("chr1.info", header=T, sep="\t", na.strings="-")
chr2 <- fread("chr2.info", header=T, sep="\t", na.strings="-")
chr3 <- fread("chr3.info", header=T, sep="\t", na.strings="-")
chr4 <- fread("chr4.info", header=T, sep="\t", na.strings="-")
chr5 <- fread("chr5.info", header=T, sep="\t", na.strings="-")
chr6 <- fread("chr6.info", header=T, sep="\t", na.strings="-")
chr7 <- fread("chr7.info", header=T, sep="\t", na.strings="-")
chr8 <- fread("chr8.info", header=T, sep="\t", na.strings="-")
chr9 <- fread("chr9.info", header=T, sep="\t", na.strings="-")
chr10 <- fread("chr10.info", header=T, sep="\t", na.strings="-")
chr11 <- fread("chr11.info", header=T, sep="\t", na.strings="-")
chr12 <- fread("chr12.info", header=T, sep="\t", na.strings="-")
chr13 <- fread("chr13.info", header=T, sep="\t", na.strings="-")
chr14 <- fread("chr14.info", header=T, sep="\t", na.strings="-")
chr15 <- fread("chr15.info", header=T, sep="\t", na.strings="-")
chr16 <- fread("chr16.info", header=T, sep="\t", na.strings="-")
chr17 <- fread("chr17.info", header=T, sep="\t", na.strings="-")
chr18 <- fread("chr18.info", header=T, sep="\t", na.strings="-")
chr19 <- fread("chr19.info", header=T, sep="\t", na.strings="-")
chr20 <- fread("chr20.info", header=T, sep="\t", na.strings="-")
chr21 <- fread("chr21.info", header=T, sep="\t", na.strings="-")
chr22 <- fread("chr22.info", header=T, sep="\t", na.strings="-")

# Bind all together
allchrm <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22)
nrow(allchrm) #Total number imputed = n
table(allchrm$Genotyped) 
# Genotyped   Imputed 
#    n  n
mean(allchrm$Rsq) #Mean of Rsq = 0.x
fivenum(allchrm$Rsq) #Five number summary of Rsq = 0.00000 0.xx 0.xx 0.xx 1.00000

## Create filter to select markers that have r2 > .3, .5, .7
allchrm_r2.3 <- allchrm[allchrm$Rsq > .3,]
nrow(allchrm_r2.3) #Total number that pass = n variants
fwrite(allchrm_r2.3, "allchrm_keep_rsq3_data.csv", row.names=F, col.names = F, quote = F)

allchrm_r2.5 <- allchrm[allchrm$Rsq > .5,]
nrow(allchrm_r2.5) #Check to get total number that pass = n
fwrite(allchrm_r2.5, "allchrm_keep_rsq5_data.csv", row.names=F, col.names = F, quote = F)

allchrm_r2.7 <- allchrm[allchrm$Rsq > .7,]
nrow(allchrm_r2.7) #Check to get total number that pass = n
fwrite(allchrm_r2.7, "allchrm_keep_rsq7_data.csv", row.names=F, col.names = F, quote = F)

allchrm_r2.9 <- allchrm[allchrm$Rsq > .9,]
nrow(allchrm_r2.9) #Check to get total number that pass = n
fwrite(allchrm_r2.9, "allchrm_keep_rsq9_data.csv", row.names=F, col.names = F, quote = F)


## Create a variable based on chrm:pos to match up to the list of HRC markers in the next step
allchrm_r2.3_chrmpos <- separate(allchrm_r2.3, SNP, c("chr", "pos"), sep = ":", extra='drop')
chrmpos_r2.3 <- paste(allchrm_r2.3_chrmpos$chr, allchrm_r2.3_chrmpos$pos, sep=":")

allchrm_r2.5_chrmpos <- separate(allchrm_r2.5, SNP, c("chr", "pos"), sep = ":", extra='drop')
chrmpos_r2.5 <- paste(allchrm_r2.5_chrmpos$chr, allchrm_r2.5_chrmpos$pos, sep=":")

allchrm_r2.7_chrmpos <- separate(allchrm_r2.7, SNP, c("chr", "pos"), sep = ":", extra='drop')
chrmpos_r2.7 <- paste(allchrm_r2.7_chrmpos$chr, allchrm_r2.7_chrmpos$pos, sep=":")

allchrm_r2.9_chrmpos <- separate(allchrm_r2.9, SNP, c("chr", "pos"), sep = ":", extra='drop')
chrmpos_r2.9 <- paste(allchrm_r2.9_chrmpos$chr, allchrm_r2.9_chrmpos$pos, sep=":")

# create update-name file for .bim file so it matches to HRC (post-imp fmt: chr:pos:a1:a2; HRC fmt: chr:pos)
allchrm$chrpos <- substr(allchrm$SNP,1,nchar(allchrm$SNP)-4) #removes the last four characters from the snp column (:a1:a2)
head(allchrm)
update_names <- allchrm[,c("SNP","chrpos")]
fwrite(update_names, "update_names.txt", sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)


## Create filter to include only the well-imputed markers that are included on the list of HRC Biallelic variants
list_HRC <- fread("/projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC_FILES/HRC_reference_biallelic_snps_bim_format_CHRPOS.txt", header=F, sep=" ") #Read list of HRC markers 
#note, the file above is the HRC file with chr:pos variant ids
#for rsid variant id format: /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC_FILES/HRC_reference_biallelic_snps_bim_format.txt
# Read 33,539,712 rows and 3 (of 3) columns from 0.258 GB file in 00:00:12

keepers_rsq3 <- list_HRC[list_HRC$V2 %in% chrmpos_r2.3,] #Creates list of markers in the list of 12mill that pass r2 threshold
nrow(keepers_rsq3) #Total number that pass = n
fwrite(keepers_rsq3[,2], "keeper_markers_rsq3_data.txt", row.names=F, col.names = F, quote = F, sep = " ")

keepers_rsq5 <- list_HRC[list_HRC$V2 %in% chrmpos_r2.5,] #Creates list of markers in the list of 12mill that pass r2 threshold
nrow(keepers_rsq5) #Total number that pass = n
fwrite(keepers_rsq5[,2], "keeper_markers_rsq5_data.txt", row.names=F, col.names = F, quote = F, sep = " ")

keepers_rsq7 <- list_HRC[list_HRC$V2 %in% chrmpos_r2.7,] #Creates list of markers in the list of 12mill that pass r2 threshold
nrow(keepers_rsq7) #Total number that pass = n
fwrite(keepers_rsq7[,2], "keeper_markers_rsq7_data.txt", row.names=F, col.names = F, quote = F, sep = " ")

keepers_rsq9 <- list_HRC[list_HRC$V2 %in% chrmpos_r2.9,] #Creates list of markers in the list of 12mill that pass r2 threshold
nrow(keepers_rsq9) #Total number that pass = 8051657
fwrite(keepers_rsq9[,2], "keeper_markers_rsq9_data.txt", row.names=F, col.names = F, quote = F, sep = " ")

## Summarize filtered SNPs:
all_imputed <- as.data.frame(table(allchrm$Genotyped))
names(all_imputed) <- c("QC", "SNPs")

rsq <- c(length(chrmpos_r2.3), length(chrmpos_r2.5), length(chrmpos_r2.7), length(chrmpos_r2.9))
rsq_label <- c("Imputation Quality: R2 > .3", "Imputation Quality: R2 > .5", "Imputation Quality: R2 > .7", "Imputation Quality: R2 > .9")
rsq2_table <- as.data.frame(cbind(rsq_label, rsq))
names(rsq2_table) <- c("QC", "SNPs")

biallelic <- c(nrow(keepers_rsq3), nrow(keepers_rsq5), nrow(keepers_rsq7), nrow(keepers_rsq9))
biallelic_label <- c("Biallelic: R2 > .3", "Biallelic: R2 > .5", "Biallelic: R2 > .7", "Biallelic: R2 > .9")
biallelic_table <- as.data.frame(cbind(biallelic_label, biallelic))
names(biallelic_table) <- c("QC", "SNPs")

imputed_table <- rbind(all_imputed, rsq2_table, biallelic_table)
write.csv(imputed_table, "data_imputed_ancestry_summary_table.csv", row.names = F, quote = T)


##### Step 4: Merge Imputed Files #####
##### Step 4: Merge Imputed Files #####

## Convert .vcf files to plink 
## Note: Because of the ID structure differences in .vcf and plink files, add "--const-fid" to this step 
## This forces zeros to be put in for FID and we will go back correct IDs in a later step
chrm <- 1:22
x3 <- NULL
for(i in 1:22) {
  x3 <- rbind(x3,paste("plink --vcf chr",chrm[i],".dose.vcf.gz --const-fid --make-bed --out data_ancestry_imputed_chr",chrm[i],sep=""))
}
write.table(x3,file="chrmerge_data_ancestry.sh",row.names=FALSE,col.names=FALSE,quote=FALSE)
system("chmod +x chrmerge_data_ancestry.sh")
system("nohup ./chrmerge_data_ancestry.sh &")


## Create list of chrm files to merge
chrm <- 1:22
x4 <- NULL
for(i in 1:22) {
  x4 <- rbind(x4,paste("data_ancestry_imputed_chr",chrm[i],sep=""))
}
write.table(x4,file="chrmerge_list.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

## Merge all 22 chrm files together
system("nohup plink --bfile data_ancestry_imputed_chr1 --merge-list chrmerge_list.txt --make-bed --out data_ancestry_imputed_allchrm --threads 20 &")
# k variants and n people pass filters and QC.

## Update names to match HRC variant ids
system("plink --bfile data_ancestry_imputed_allchrm --update-name update_names.txt --make-bed --out allchrm_data_ancestry_chrpos")
# k variants and n people pass filters and QC.

##### Step 5: Update participant IDs to match original IDs #####
##### Step 5: Update participant IDs to match original IDs #####

## Update participant IDs to match original .fam file
## Read in original .fam file
original_fam <- fread("filepath-to-original-data/data.fam")
head(original_fam)

table(original_fam$V1) ## Check if unique FAM ID's exist 

## Read in imputed .fam file
imp_fam <- fread("filepath-to-imputed-data/allchrm_data_ancestry_chrpos.fam", header=F, sep=" ")
head(imp_fam)

## Need to split V2 into FID/IID, replace _ with , and split into new variables
imp_fam$newV2 <- gsub("_", ",", imp_fam$V2)
imp_fam_new <- imp_fam %>% separate(col = newV2, into = c("FID","IID"), sep = ",")

## Write new .fam files
write.table(imp_fam_new[,c(7:8,3:6)], "imp_fam_IDadjust.fam", row.names = F, col.names = F, quote = F) 

## Update the IDs in the imputed file
system("plink --bfile allchrm_data_ancestry_chrpos --fam imp_fam_IDadjust.fam --make-bed --out allchrm_data_ancestry_chrpos_ids")
# k variants and n people pass filters and QC.


##### Step 6: Screen Markers #####
##### Step 6: Screen Markers #####

# Screen for the biallelic markers that were well imputed (rsq > .3 and present in HRC file)
system("plink --bfile allchrm_data_ancestry_chrpos_ids --extract keeper_markers_rsq3_data_ancestry.txt --make-bed --out data_ancestry_imputed_rsq3")
# k variants and n people pass filters and QC.

#screen for the biallelic markers that were well imputed (rsq > .7 and present in HRC file)
system("plink --bfile allchrm_data_ancestry_chrpos_ids --extract keeper_markers_rsq7_data_ancestry.txt --make-bed --out data_ancestry_imputed_rsq7")
# k variants and n people pass filters and QC.


## Check for duplicate markers in the imputed data
bim <- fread("data_ancestry_imputed_rsq3.bim")
head(bim)
table(duplicated(bim$V2)) ## chr:pos
# FALSE 
# k

## Create list of unique variants to keep (this removes BOTH occurrances of a SNP)
keepers <- bim[!(duplicated(bim$V2) | duplicated(bim$V2, fromLast = TRUE)), ]
dim(keepers)
#  k        6

dups <- bim[(duplicated(bim$V2) | duplicated(bim$V2, fromLast = TRUE)), ]
dim(dups)
# 0 6

## Write out list of makers to screen because they are duplicated
fwrite(keepers[,2], "unique_markers_rsq3.txt", row.names=F, col.names = F, quote = F, sep = " ")
fwrite(dups, "dup_markers_rsq3.txt", row.names=F, col.names = F, quote = F, sep = " ")

## Screen duplicates out of imputed plink file
system("plink --bfile data_ancestry_imputed_rsq3 --extract unique_markers_rsq3.txt  --make-bed --out data_ancestry_imputed_rsq3_nodup > data_ancestry_imputed_rsq3_nodup.txt")
# k variants and n people pass filters and QC.

#### FINAL FILE = data_ancestry_imputed_rsq3_nodup #### 
#### FINAL FILE = data_ancestry_imputed_rsq3_nodup #### ==> NEED TO CONVERT TO FORMAL NAMING STRUCTURE WHEN FINALIZED!
#### FINAL FILE = data_ancestry_imputed_rsq3_nodup #### 
system("plink --bfile data_ancestry_imputed_rsq3  --make-bed --out FINAL_NAME_TBD")