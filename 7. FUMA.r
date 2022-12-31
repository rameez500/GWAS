assoc <- fread("../MLMA/post_imp_factor_res_age_sex_PC14_cov_cohort_ethnic_mlma.mlma")
assoc <- assoc[order(assoc$p),]
head(assoc,30); dim(assoc)

#(1) Convert your output to chi-squared values

# For p-values, calculate chi-squared statistic
chisq <- qchisq(1-assoc$p,1)

#(2) Calculate lambda gc (λgc)
median(chisq)/qchisq(0.5,1)
# 0.9882559

# If lambda is less than or equal to 1, no adjustment is necessary
# If the λ value is greater than 1, then this may be evidence for some systematic bias that needs to be corrected in your analysis


fuma <- assoc[,c("Chr","bp","p")]
colnames(fuma)[1:3] <- c("CHR","BP","P")
head(fuma)

fwrite(fuma,"fuma_load.txt",sep = "\t")


# Result from FUMA 
setwd("FUMA_job220960/")
annot <- fread("annot.txt")
annot <- annot[,c(1:3)]
head(annot); dim(annot) 


annov_stat <-fread("annov.stats.txt")
head(annov_stat); dim(annov_stat)


annov <- fread("annov.txt")
annov$chr_pos <- paste(annov$chr,annov$pos,sep = ":")
head(annov); dim(annov)

genes <- fread("genes.txt")
head(genes); dim(genes)

genomicriskloci <- fread("GenomicRiskLoci.txt")
head(genomicriskloci); dim(genomicriskloci)

gwascatolog <- fread("gwascatalog.txt")
head(gwascatolog); dim(gwascatolog)

indsigSNP <- fread("IndSigSNPs.txt")
head(indsigSNP); dim(indsigSNP)

id <- fread("ld.txt")
head(id); dim(id)

leadSNP <- fread("leadSNPs.txt")
head(leadSNP)


magma_genes <-  fread("magma.genes.out")
head(magma_genes); dim(magma_genes)

