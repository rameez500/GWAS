####################  H-magma ##################################
####################  H-magma ##################################
####################  H-magma ##################################


assoc <- fread("../MLMA/post_imp_factor_res_age_sex_PC14_cov_cohort_ethnic_mlma.mlma")
assoc <- assoc[order(assoc$p),]
head(assoc,30); dim(assoc)

bim <- fread("rs_1k.txt")
bim <- bim[!duplicated(bim$V7),]

bim_assoc <- merge(assoc,bim,by.x = "SNP", by.y = "V7", sort = F)

bim_assoc_sub <- bim_assoc[,c("V2","p")]
colnames(bim_assoc_sub)[1:2] <- c("SNP","P")
head(bim_assoc_sub); dim(bim_assoc_sub)

fwrite(bim_assoc_sub[,c("SNP","P")],"sumStat_rs.txt",sep = "\t",col.names = TRUE)

# Gene analysis on SNP p-value data
# magma default
/home/rasyed2/Tool/magma --bfile /projects/bga_lab/DATA_REPOSITORIES/post_imp_WITHDRAWAL/UM_IMP_RESULTS/MERGED_CLEAN/SmokeScreen/Association_V3/1000_mix/g1000Mix.QC2 \
--pval sumStat_rs.txt N=8009 \
--gene-annot /projects/bga_lab/DATA_REPOSITORIES/post_imp_WITHDRAWAL/UM_IMP_RESULTS/MERGED_CLEAN/SmokeScreen/Association_V3/H_magma/MAGMAdefault.genes.annot \
--out gene_based_H_MAGMAdefault                                

# adult 
/home/rasyed2/Tool/magma --bfile /projects/bga_lab/DATA_REPOSITORIES/post_imp_WITHDRAWAL/UM_IMP_RESULTS/MERGED_CLEAN/SmokeScreen/Association_V3/1000_mix/g1000Mix.QC2 \
--pval sumStat_rs.txt N=8009 \
 --gene-annot /projects/bga_lab/DATA_REPOSITORIES/post_imp_WITHDRAWAL/UM_IMP_RESULTS/MERGED_CLEAN/SmokeScreen/Association_V3/H_magma/Adult_brain.genes.annot \
 --out gene_based_H_Adult_brain                                

# fetal   
/home/rasyed2/Tool/magma --bfile /projects/bga_lab/DATA_REPOSITORIES/post_imp_WITHDRAWAL/UM_IMP_RESULTS/MERGED_CLEAN/SmokeScreen/Association_V3/1000_mix/g1000Mix.QC2 \
--pval sumStat_rs.txt N=8009 \
--gene-annot /projects/bga_lab/DATA_REPOSITORIES/post_imp_WITHDRAWAL/UM_IMP_RESULTS/MERGED_CLEAN_V2/gene_based/H_magma/Fetal_brain.genes.annot \
--out gene_based_H_Fetal_brain                                


# neuronal  
/home/rasyed2/Tool/magma --bfile /projects/bga_lab/DATA_REPOSITORIES/post_imp_WITHDRAWAL/UM_IMP_RESULTS/MERGED_CLEAN/SmokeScreen/Association_V3/1000_mix/g1000Mix.QC2 \
--pval sumStat_rs.txt N=8009 \
--gene-annot /projects/bga_lab/DATA_REPOSITORIES/post_imp_WITHDRAWAL/UM_IMP_RESULTS/MERGED_CLEAN_V2/gene_based/H_magma/' iPSC_derived_neuro.genes.annot' \
--out gene_based_H_neuro                                


# astrocyte  
/home/rasyed2/Tool/magma --bfile /projects/bga_lab/DATA_REPOSITORIES/post_imp_WITHDRAWAL/UM_IMP_RESULTS/MERGED_CLEAN/SmokeScreen/Association_V3/1000_mix/g1000Mix.QC2 \
--pval sumStat_rs.txt N=8009 \
--gene-annot /projects/bga_lab/DATA_REPOSITORIES/post_imp_WITHDRAWAL/UM_IMP_RESULTS/MERGED_CLEAN_V2/gene_based/H_magma/iPSC_derived_astro.genes.annot  \
--out gene_based_H_astro                                


adult <- fread("gene_based_H_Adult_brain.genes.out")
adult <- adult[,c("GENE","ZSTAT","P")]
colnames(adult)[2:3] <- c("Z_adult","P_adult")
head(adult); dim(adult)

fetal <- fread("gene_based_H_Fetal_brain.genes.out")
fetal <- fetal[,c("GENE","ZSTAT","P")]
colnames(fetal)[2:3] <- c("Z_fetal","P_fetal")
head(fetal); dim(fetal)


astrocyte <- fread("gene_based_H_astro.genes.out")
astrocyte <- astrocyte[,c("GENE","ZSTAT","P")]
colnames(astrocyte)[2:3] <- c("Z_astro","P_astro")
head(astrocyte); dim(astrocyte)


neuronal <- fread("gene_based_H_neuro.genes.out")
neuronal <- neuronal[,c("GENE","ZSTAT","P")]
colnames(neuronal)[2:3] <- c("Z_neuro","P_neuro")
head(neuronal); dim(neuronal)


adult_fetal <- merge(adult,fetal, by = "GENE", sort = FALSE)
head(adult_fetal); dim(adult_fetal)

adult_fetal_astro <- merge(adult_fetal,astrocyte, by = "GENE", sort = FALSE)
head(adult_fetal_astro); dim(adult_fetal_astro)

adult_fetal_astro_neuro <- merge(adult_fetal_astro,neuronal, by = "GENE", sort = FALSE)
head(adult_fetal_astro_neuro); dim(adult_fetal_astro_neuro)


#### option 1 ####
library(biomaRt)
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
gene_37 <-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version","hgnc_symbol",'chromosome_name','start_position','end_position'), mart = grch37)
gene_37 <- gene_37[gene_37$chromosome_name == 1 | gene_37$chromosome_name == 2 | gene_37$chromosome_name == 3 | gene_37$chromosome_name == 4 |
				  gene_37$chromosome_name == 5 | gene_37$chromosome_name == 6 | gene_37$chromosome_name == 7 | gene_37$chromosome_name == 8  |
				  gene_37$chromosome_name == 9 | gene_37$chromosome_name == 10 | gene_37$chromosome_name == 11 | gene_37$chromosome_name == 12  |
				  gene_37$chromosome_name == 13 | gene_37$chromosome_name == 14 | gene_37$chromosome_name == 15 | gene_37$chromosome_name == 16  |
				  gene_37$chromosome_name == 17 | gene_37$chromosome_name == 18 | gene_37$chromosome_name == 19 | gene_37$chromosome_name == 20  |
				  gene_37$chromosome_name == 21 | gene_37$chromosome_name == 22 ,]

head(gene_37); dim(gene_37)
table(gene_37$chromosome_name)

hmagma <- merge(adult_fetal_astro_neuro,gene_37,by.x = "GENE", by.y = "ensembl_gene_id", sort = F, all.x = TRUE)
hmagma <- hmagma[,c(11,1:9)]
hmagma$hgnc_symbol <- as.character(hmagma$hgnc_symbol)
hmagma <- as.data.frame(hmagma)
head(hmagma); dim(hmagma); dim(hmagma[!duplicated(hmagma$GENE),])

sum(is.na(hmagma$hgnc_symbol))

head(hmagma[!is.na(hmagma$hgnc_symbol),])

table((hmagma$hgnc_symbol == ""))


#### option 2 ####
# https://www.biotools.fr/human/ensembl_symbol_converter
fwrite(as.data.frame(hmagma[,2]),"ENS_gene.txt", col.names = F)

convert_eng_gene <- fread("convert_eng_gene.txt",header = FALSE)
head(convert_eng_gene); dim(convert_eng_gene)

hmagma <- merge(hmagma,convert_eng_gene,by.x = "GENE", by.y = "V1", sort = FALSE)
hmagma <- hmagma[!duplicated(hmagma$GENE),]
head(hmagma); dim(hmagma)

table(hmagma$hgnc_symbol == hmagma$V2)

hmagma <- hmagma[,c(1,11,3:10)]
colnames(hmagma)[2] <- "hgnc_symbol"
head(hmagma); dim(hmagma)



library(matrixStats)
 
df1 <- rowMins(as.matrix((hmagma[,c(4,6,8,10)]) ))
hmagma$P_Minimum_of_adult_fetal_astro_neuro <- df1

#  FDR
Padj_of_P_Minimum <- p.adjust(hmagma$P_Minimum_of_adult_fetal_astro_neuro, method = "fdr", n = length(hmagma$P_Minimum_of_adult_fetal_astro_neuro))
hmagma <- cbind(hmagma,Padj_of_P_Minimum)
hmagma <- as.data.frame(hmagma)
hmagma <- hmagma[order(hmagma$P_Minimum_of_adult_fetal_astro_neuro),]
