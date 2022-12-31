pheno_mod <- lm(Factore_Score ~ Gender.x +  Age.x +  PC1 + PC2 + PC3 + PC4 ,
 data = Merge_geno_pheno_fam_covar)
summary(pheno_mod)


new_standard_res_PC1_4_study_cov <- rstandard(pheno_mod)
cohort_post_imp_geno_pheno <- cbind(Merge_geno_pheno_fam_covar,new_standard_res_PC1_4_study_cov)
cohort_post_imp_geno_pheno <- as.data.frame(cohort_post_imp_geno_pheno)
head(cohort_post_imp_geno_pheno); dim(cohort_post_imp_geno_pheno)

fwrite(cohort_post_imp_geno_pheno,"MLMA/post_imp_mix_pheno.csv")

write.table(cohort_post_imp_geno_pheno[,c(2,1,33)],
"MLMA/factor_pheno_res_PC14_age_sex_cov.txt", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t", na = "NA")


