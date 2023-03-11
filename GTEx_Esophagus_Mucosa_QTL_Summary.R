#! /usr/lib/R/bin/Rscript
####################################################################################################
##README:
##This script aims to
#1) Summary the eQTL mapping results
#2) compare the eQTL results with different bins
####################################################################################################
library(data.table)
##step I: Load the ieQTL_(top QTL) results (in R)
bin1 = as.data.frame(fread("../GTEx_Esophagus_Mucosa_bin1.cis_qtl_top_assoc.txt.gz"))
bin2 = as.data.frame(fread("../GTEx_Esophagus_Mucosa_bin2.cis_qtl_top_assoc.txt.gz"))
bin3 = as.data.frame(fread("../GTEx_Esophagus_Mucosa_bin3.cis_qtl_top_assoc.txt.gz"))
bin4 = as.data.frame(fread("../GTEx_Esophagus_Mucosa_bin4.cis_qtl_top_assoc.txt.gz"))
#top_eQTL of each gene (across bins)
p_adj = cbind(bin4$pval_adj_bh,bin3$pval_adj_bh,bin2$pval_adj_bh,bin1$pval_adj_bh)
p_gi = cbind(bin4$pval_gi,bin3$pval_gi,bin2$pval_gi,bin1$pval_gi)
pval_emt = cbind(bin4$pval_emt,bin3$pval_emt,bin2$pval_emt,bin1$pval_emt)
variant = cbind(bin4$variant_id,bin3$variant_id,bin2$variant_id,bin1$variant_id)
Tss = cbind(bin4$tss_distance,bin3$tss_distance,bin2$tss_distance,bin1$tss_distance)
maf = cbind(bin4$af,bin3$af,bin2$af,bin1$af)

topBin = apply(p_adj,1,which.min)
top_variant = sapply(1:length(topBin),function(i){variant[i,topBin[i]]})
top_p_gi = sapply(1:length(topBin),function(i){p_gi[i,topBin[i]]})
top_Tss = sapply(1:length(topBin),function(i){Tss[i,topBin[i]]})
top_maf = sapply(1:length(topBin),function(i){maf[i,topBin[i]]})
top_p = sapply(1:length(topBin),function(i){p_adj[i,topBin[i]]})
top_pval_emt = sapply(1:length(topBin),function(i){pval_emt[i,topBin[i]]})
commonGene = bin4$phenotype_id

out = data.frame(variant_id=top_variant,gene_id=commonGene,maf=top_maf,tss_distance=top_Tss,pval_emt=top_pval_emt,pval_adj_bh=top_p)
write.table(out,"../Esophagus_Mucosa.cell_state.ieQTL.eigenMT.annotated.txt",
	col.names=T,sep='\t',quote=F)
#permute the p-adj, prepare for following analyses
out_random = data.frame(variant_id=top_variant,gene_id=commonGene,maf=top_maf,tss_distance=top_Tss,pval_emt=top_pval_emt,pval_adj_bh=sample(top_p))
write.table(out_random,"../Esophagus_Mucosa.cell_state.ieQTL.eigenMT.random.annotated.txt",
	col.names=T,sep='\t',quote=F)

saveRDS(list(bin1,bin2,bin3,bin4),"../top_ieQTL_Esophagus_Mucosa.rds")


##step II: compare the eQTL results with different bins (in R)
#using the FDR threshold < 0.05
bin1 = bin1[bin1$pval_adj_bh<0.05,]
bin2 = bin2[bin2$pval_adj_bh<0.05,]
bin3 = bin3[bin3$pval_adj_bh<0.05,]
bin4 = bin4[bin4$pval_adj_bh<0.05,]
#eGene
eGene = unique(c(bin1$phenotype_id,bin2$phenotype_id,bin3$phenotype_id,bin4$phenotype_id))
#Load the data of 10bins
sensitivity = sapply(1:10,function(i){
	path = paste0("../GTEx_Esophagus_Mucosa",paste0(paste0("_bin",i),".sensitivity.cis_qtl_top_assoc.txt.gz"))
	temp = as.data.frame(fread(path))
	list(temp)
})
#Compare the p-value of eGene with other Genes in 10 bins
sensitivity_p = do.call(cbind,sapply(sensitivity,function(i){list(i$pval_gi)}))
topBin_sensitivity = apply(sensitivity_p,1,which.min)
sensitivity_p = sapply(1:length(topBin_sensitivity),function(i){sensitivity_p[i,topBin_sensitivity[i]]})
geneType = rep('eGene',length(sensitivity_p))
geneType[!sensitivity[[1]]$phenotype_id %in% eGene] = 'otherGene'
sensitivity_analysis = data.frame(top_p_gi,sensitivity_p,geneType)
saveRDS(sensitivity_analysis,'../sensitivity_analysis.rds')