#  FDR
p_adjust_fdr <- p.adjust(magma_genes$P, method = "fdr", n = length(magma_genes$P))
# p_adjust_fdr <- p.adjust(magma_genes$P, method = "BH", n = length(magma_genes$P))
magma_genes <- cbind(magma_genes,p_adjust_fdr)
magma_genes <- as.data.frame(magma_genes)
magma_genes <- magma_genes[order(magma_genes$P),]
magma_genes <- magma_genes[,c(10,2:5,8,9,11)]
head(magma_genes) ; dim(magma_genes)
