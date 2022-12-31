################ 	Manhattan Plot  	####################
################ 	Manhattan Plot  	####################
################ 	Manhattan Plot  	####################

# Intergenic

mlma_fuma_annovar_noDup_intergenic <- sig_assoc_annot[sig_assoc_annot$annot == "intergenic" | 
                                                              sig_assoc_annot$annot == "ncRNA_intronic" |
                                                              sig_assoc_annot$annot == "ncRNA_exonic",]

mlma_fuma_annovar_noDup_intergenic; dim(mlma_fuma_annovar_noDup_intergenic)


intergenic_snp <- mlma_fuma_annovar_noDup_intergenic[,c("SNP","annot")]
colnames(intergenic_snp)[1:2] <- c("SNP","Genes")
intergenic_snp



# Intronic/exonic

mlma_fuma_annovar_noDup_intronic <- sig_assoc_annot[sig_assoc_annot$annot == "intronic" 
                                                              | sig_assoc_annot$annot == "exonic",]

mlma_fuma_annovar_noDup_intronic; dim(mlma_fuma_annovar_noDup_intronic)

intron_snp <- mlma_fuma_annovar_noDup_intronic[,c("SNP","symbol")]
colnames(intron_snp)[1:2] <- c("SNP","Genes")
intron_snp



intro_intergenic <- rbind(intergenic_snp,intron_snp)
intro_intergenic

SNPs <- as.factor(intro_intergenic$SNP)
SNPs


genes <- intro_intergenic$Genes
genes


x <- assoc[,c("SNP","Chr","bp","p")]
colnames(x) <- c("SNP","Chromosome","Position","p")
head(x); dim(x)


CMplot(x, plot.type="m", col=c("blue","grey30"), LOG10=TRUE, ylim=c(0,28), 
       threshold=c(5e-8,1e-5),cex = 0.5, highlight.cex = 0.5,
        threshold.lty=c(1,2), threshold.lwd=c(1,2), threshold.col=c("red","black"), amplify=FALSE,highlight=SNPs,
         highlight.text=genes, highlight.col= NULL, file="jpg",memo="coga_mix",dpi=600,file.output=TRUE,verbose=TRUE,width=12,height=12)
		 