########## Combime assoc and FUMA Table ####################
########## Combime assoc and FUMA Table ####################
########## Combime assoc and FUMA Table ####################

assoc_annov <- merge(assoc,annov,by.x = "SNP", by.y = "chr_pos", sort = FALSE)
head(assoc_annov); dim(assoc_annov)


assoc_annov_annot <- merge(assoc_annov,annot,by = "uniqID", sort = FALSE)
assoc_annov_annot <- assoc_annov_annot[,c("uniqID", "Chr","bp", "Freq", "p", "b","A1", "symbol","annot", "exonic_func", "CADD", "RDB","SNP")]
head(assoc_annov_annot); dim(assoc_annov_annot)


fwrite(assoc_annov_annot,"../mlma_fuma_annovar.csv", sep = ",")
# merge with 1000G to get rsid

assoc_annov_annot_rsid <- merge(assoc_annov_annot,bim,by.x = "SNP",by.y = "V7",all.x = TRUE)
assoc_annov_annot_rsid <- assoc_annov_annot_rsid[,c("V2","Chr","bp","A1","Freq","p"
											,"b","symbol","annot","CADD","RDB","SNP")]
assoc_annov_annot_rsid <- assoc_annov_annot_rsid[order(assoc_annov_annot_rsid$p),]											
head(assoc_annov_annot_rsid); dim(assoc_annov_annot_rsid)


fwrite(assoc_annov_annot_rsid,"../mlma_fuma_annovar_V2.csv", sep = ",")

#############  genome wide association threshold (5 X 10^-8)   #######################
#############  genome wide association threshold (5 X 10^-8)   #######################
#############  genome wide association threshold (5 X 10^-8)   #######################

assoc_annov_annot_gwas_thres <- assoc_annov_annot[assoc_annov_annot$p < 5e-8,]
head(assoc_annov_annot_gwas_thres); dim(assoc_annov_annot_gwas_thres)



#############  suggestive  threshold (1 X 10^-5)   #######################
#############  suggestive  threshold (1 X 10^-5)   #######################
#############  suggestive  threshold (1 X 10^-5)   #######################


assoc_annov_annot_sugg_thres <- assoc_annov_annot[assoc_annov_annot$p >= 5e-8 & assoc_annov_annot$p < 1e-5,]
head(assoc_annov_annot_sugg_thres); dim(assoc_annov_annot_sugg_thres)



sig_assoc_annot <- assoc_annov_annot_rsid[assoc_annov_annot_rsid$p < 1e-5,]
sig_assoc_annot <- sig_assoc_annot[!duplicated(sig_assoc_annot$SNP),]
sig_assoc_annot <- sig_assoc_annot[!duplicated(sig_assoc_annot$symbol),]
head(sig_assoc_annot); dim(sig_assoc_annot)
fwrite(sig_assoc_annot,"../mlma_fuma_annovar_V3.csv", sep = ",")


assoc_annov_annot_rsid[(assoc_annov_annot_rsid$V2 %in% sig_assoc_annot$V2[1:2] ),]