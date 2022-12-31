#########################  MetaXcan ##############################
#########################  MetaXcan ##############################
#########################  MetaXcan ##############################

# update the SNP id for tissue expression
bim_assoc_update <- bim_assoc[,c(2,11,3:9)]
colnames(bim_assoc_update)[2] <- "SNP"
head(bim_assoc_update); dim(bim_assoc_update)
fwrite(bim_assoc_update,"../Multixcan/gwas/AWS_GWAS_summary.txt", sep = "\t")

# Predict gene expression using a single tissue


### AMYGDALA

/projects/bga_lab/SOFTWARE/Tool/MetaXcan/software/./SPrediXcan.py \
--model_db_path en_Brain_Amygdala.db \
--covariance en_Brain_Amygdala.txt.gz \
--gwas_folder gwas/ \
--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p \
--output_file brain_en_Brain_Amygdala.csv

#
# ACC

/projects/bga_lab/SOFTWARE/Tool/MetaXcan/software/./SPrediXcan.py \
--model_db_path en_Brain_Anterior_cingulate_cortex_BA24.db \
--covariance en_Brain_Anterior_cingulate_cortex_BA24.txt.gz \
--gwas_folder gwas/ \
--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p \
--output_file brain_en_Brain_Anterior_cingulate_cortex_BA24.csv

## Caudate

/projects/bga_lab/SOFTWARE/Tool/MetaXcan/software/./SPrediXcan.py \
--model_db_path en_Brain_Caudate_basal_ganglia.db --covariance en_Brain_Caudate_basal_ganglia.txt.gz \
--gwas_folder gwas/ \
--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p \
--output_file brain_en_Brain_Caudate_basal_ganglia.csv

# Cerebellar_Hemisphere

/projects/bga_lab/SOFTWARE/Tool/MetaXcan/software/./SPrediXcan.py \
--model_db_path en_Brain_Cerebellar_Hemisphere.db --covariance en_Brain_Cerebellar_Hemisphere.txt.gz \
--gwas_folder gwas/ \
--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p \
--output_file brain_en_Brain_Cerebellar_Hemisphere.csv


## Cerebellum

/projects/bga_lab/SOFTWARE/Tool/MetaXcan/software/./SPrediXcan.py \
--model_db_path en_Brain_Cerebellum.db \
--covariance en_Brain_Cerebellum.txt.gz \
--gwas_folder gwas/ \
--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p \
--output_file brain_en_Brain_Cerebellum.csv

## Cortex 

/projects/bga_lab/SOFTWARE/Tool/MetaXcan/software/./SPrediXcan.py \
--model_db_path en_Brain_Cortex.db \
--covariance en_Brain_Cortex.txt.gz \
--gwas_folder gwas/ \
--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p \
--output_file brain_en_Brain_Cortex.csv

## PFC - BA9

/projects/bga_lab/SOFTWARE/Tool/MetaXcan/software/./SPrediXcan.py \
--model_db_path en_Brain_Frontal_Cortex_BA9.db \
--covariance en_Brain_Frontal_Cortex_BA9.txt.gz \
--gwas_folder gwas/ \
--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p \
--output_file brain_en_Brain_Frontal_Cortex_BA9.csv

## Hippocampus

/projects/bga_lab/SOFTWARE/Tool/MetaXcan/software/./SPrediXcan.py \
--model_db_path en_Brain_Hippocampus.db \
--covariance en_Brain_Hippocampus.txt.gz \
--gwas_folder gwas/ \
--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p \
--output_file brain_en_Brain_Hippocampus.csv

## Hypothalamus

/projects/bga_lab/SOFTWARE/Tool/MetaXcan/software/./SPrediXcan.py \
--model_db_path en_Brain_Hypothalamus.db \
--covariance en_Brain_Hypothalamus.txt.gz \
--gwas_folder gwas/ \
--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p \
--output_file brain_en_Brain_Hypothalamus.csv

## Nucleus Accumbens

/projects/bga_lab/SOFTWARE/Tool/MetaXcan/software/./SPrediXcan.py \
--model_db_path en_Brain_Nucleus_accumbens_basal_ganglia.db \
--covariance en_Brain_Nucleus_accumbens_basal_ganglia.txt.gz \
--gwas_folder gwas/ \
--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p \
--output_file brain_en_Brain_Nucleus_accumbens_basal_ganglia.csv

## Putamen

/projects/bga_lab/SOFTWARE/Tool/MetaXcan/software/./SPrediXcan.py \
--model_db_path en_Brain_Putamen_basal_ganglia.db \
--covariance en_Brain_Putamen_basal_ganglia.txt.gz \
--gwas_folder gwas/ \
--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p \
--output_file brain_en_Brain_Putamen_basal_ganglia.csv

## Spinal Cord

/projects/bga_lab/SOFTWARE/Tool/MetaXcan/software/./SPrediXcan.py \
--model_db_path en_Brain_Spinal_cord_cervical_c-1.db \
--covariance en_Brain_Spinal_cord_cervical_c-1.txt.gz \
--gwas_folder gwas/ \
--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p \
--output_file brain_en_Brain_Spinal_cord_cervical_c-1.csv

## Substantia nigra

/projects/bga_lab/SOFTWARE/Tool/MetaXcan/software/./SPrediXcan.py \
--model_db_path en_Brain_Substantia_nigra.db \
--covariance en_Brain_Substantia_nigra.txt.gz \
--gwas_folder gwas/ \
--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p \
--output_file brain_en_Brain_Substantia_nigra.csv




# Predict gene expression using multiple tissues

/projects/bga_lab/SOFTWARE/Tool/MetaXcan/software/./SMulTiXcan.py \
--models_folder gene_expression/GTE_data/elastic_net_models \
--models_name_pattern "en_(.*).db" \
--snp_covariance rediXscan/gtex_v8_expression_elastic_net_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder  \
--metaxcan_filter "brain_en_(.*).csv" \
--metaxcan_file_name_parse_pattern "(.*)_en_(.*).csv" \
--gwas_file gwas/AWS_GWAS_summary.txt \
--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output BRAIN.txt