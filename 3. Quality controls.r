library(data.table)
library(dplyr)
library(tidyr)


## step 1: geno 0.1
system("plink --bfile post_imp --geno 0.1 --make-bed --out post_imp_geno")



## step 2: hwe 0.0001
system("plink --bfile post_imp_geno --hwe 0.0001 --make-bed --out post_imp_geno_hwe")


## step 3: maf 0.05
system("plink --bfile post_imp_geno_hwe --maf 0.05 --make-bed --out post_imp_geno_hwe_maf")


## step 4: mind 0.1
system("plink --bfile post_imp_geno_hwe_maf --mind 0.1 --make-bed --out post_imp_geno_hwe_maf_mind")









#----