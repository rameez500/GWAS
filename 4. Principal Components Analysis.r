## 1000 Genomes LD pruned data (50 10 .5) file (original file from SVS)
kg_bim <- fread("/projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/IMPUTE2_1KG_Reference_All_Biallelic_SNPS.bim")
dim(kg_bim) # 2 050 463       13
kg_bim$chr_pos <- paste(kg_bim$V1,kg_bim$V4,sep = ":")
head(kg_bim)
kg_bim_snps <- as.matrix(kg_bim$chr_pos)
sum(duplicated(kg_bim_snps))


## Read in your sample data .bim
post_imp_bim <- fread("post_imp_geno_hwe_maf_mind.bim")
head(post_imp_bim)
nrow(post_imp_bim) # n variants
post_imp_bim_snps <- as.matrix(post_imp_bim$V2)



table(kg_bim_snps %in% post_imp_bim_snps)
overlap <- kg_bim_snps[kg_bim_snps %in% post_imp_bim_snps]
length(overlap) # n
library(data.table)
fwrite(as.list(overlap), "post_imp_1kg_ref_snps.txt", sep = " ", quote = FALSE, col.names = FALSE)



system(
"plink --bfile post_imp_geno_hwe_maf_mind \\
--extract post_imp_1kg_ref_snps.txt --make-bed --out  post_imp_match_1kgRef")


system("flashpca --bfile post_imp_match_1kgRef --ndim 20 --suffix _post_imp.txt --numthreads 100")




# Data Visualization: Scree Plot
eigen_value <- fread("pve_post_imp.txt")
eigen_value <- eigen_value[1:10,]
seq_PC <- seq(1,10)
eigen_value$PC <- as.factor(paste("PC",seq_PC,sep = ""))
eigen_value <- as.data.frame(eigen_value)
head(eigen_value)

png("plot_pc.png")
ggplot(eigen_value, aes(x = fct_inorder(PC),y = V1)) + xlab("Principal Components") + ylab("Eigenvalues") + 
theme(
axis.title.x = element_text(size=14, face="bold"),
axis.title.y = element_text(size=14, face="bold"),
axis.text = element_text(size = (13))
) + 
  geom_point()
dev.off()