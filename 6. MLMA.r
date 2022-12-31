
gcta193 --bfile post_imp_WD_imp_all_gr_maf_mind_pheno --autosome --thread-num 50 --make-grm-bin --out post_imp_WD
gcta193 --grm post_imp_WD --make-bK 0.05 --out post_imp_WDgrm_bk --thread-num 50 

echo post_imp_WD >> post_imp_WD_NormalNames_IMP_GRMs.txt
echo post_imp_WDgrm_bk >> post_imp_WD_NormalNames_IMP_GRMs.txt



gcta193 --mlma --bfile post_imp_geno_hwe_maf_mind \
--mgrm post_imp_WD_NormalNames_IMP_GRMs.txt \
--pheno factor_pheno_res_PC14_age_sex_cov.txt \
--thread-num 80 \
--out post_imp_factor_res_age_sex_PC14_cov_cohort_ethnic_mlma


gcta193 --reml \
--mgrm post_imp_WD_NormalNames_IMP_GRMs.txt \
--reml-maxit 200 \
--pheno factor_pheno_res_PC14_age_sex_cov.txt \
--thread-num 70 \
--out post_imp_factor_res_age_sex_PC14_cov_cohort_ethnic_reml
