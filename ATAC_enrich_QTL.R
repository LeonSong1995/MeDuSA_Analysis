#! /usr/lib/R/bin/Rscript
###Permutation test of the scATAC-seq
library(data.table)
library(argparser)
library(viridisLite)
library(ggplot2)
library(parallel)
library(GenomicRanges)

p = arg_parser("Prepare SMR esd")
p = add_argument(p, "ieQTL", help="")
p = add_argument(p, "out", help="")
p = add_argument(p, "peak", help="")
p = add_argument(p, "p", help="")
p = add_argument(p, "type", help="")
argv = parse_args(p)

type = argv$type

#null: 1 is null o is eQTL
null = as.data.frame(fread("../GTEx_v8_Null_Table_49Tiss_Global.txt.gz"))
confounders = as.data.frame(readRDS("../GTEx_v8_Confounders_Table_49Tiss_Global.bined.rds"))
ieQTL = as.data.frame(fread(argv$ieQTL))
if(type=='eQTL'){
	ieQTL = ieQTL[ieQTL$pval_beta<argv$p,]
}else if(type=='ieQTL'){
	ieQTL = ieQTL[ieQTL$pval_adj_bh<argv$p,]
}
ieQTL$chr = sapply(ieQTL$variant_id,function(i){strsplit(i,split = '_')[[1]][1]})
ieQTL$pos = as.numeric(sapply(ieQTL$variant_id,function(i){strsplit(i,split = '_')[[1]][2]}))
peak = readRDS(argv$peak)
peak = peak[peak$score>0,]


###prase the ieQTL table
index = match(ieQTL[,'variant_id'],confounders[,'variant'])
ieQTL$MAF_Bin = confounders$MAF_Bin[index]
ieQTL$TSS_Bin = confounders$TSS_Bin[index]
ieQTL$LD_Bin = confounders$LD_Bin[index]

##Run permutations
null_variat = c(null$variant[null$Global==1])
null_variat = null_variat[!null_variat %in% ieQTL$variant_id]
confounders = confounders[confounders$variant %in% null_variat,]
subset = sapply(1:nrow(ieQTL),function(i){
	print(i)
	temp = ieQTL[i,]
	eligible = Reduce(intersect,list(which(temp$MAF_Bin==confounders$MAF_Bin),
									 which(temp$TSS_Bin==confounders$TSS_Bin),
									 which(temp$LD_Bin==confounders$LD_Bin)))
	list(eligible)
})
subset = subset[lapply(subset,length)!=0]

permute = sapply(1:10000,function(i){
	print(i)
	id = sapply(subset,function(i){sample(i,1)})
	pseudo = confounders[id,]
	pseudoQTL = GRanges(seqnames = paste0(paste('b',pseudo$chr,sep="'"),"'"),
	               ranges = IRanges(start = pseudo$pos-100,end = pseudo$pos+100))
	r=peak[GenomicRanges::findOverlaps(pseudoQTL,peak)@to,]$score
	mean(r)
})
permute = na.omit(permute)


csd_eQTL = GRanges(seqnames = paste0(paste('b',ieQTL$chr,sep="'"),"'"),
               ranges = IRanges(start = ieQTL$pos-100,end = ieQTL$pos+100))


observed_fold_enrichment = mean(peak[GenomicRanges::findOverlaps(csd_eQTL,peak)@to,]$score)
first_fold_enrichments = permute
median_fold_enrichment = median(first_fold_enrichments)
adjusted_fold_enrichment = round(observed_fold_enrichment/median_fold_enrichment,4)


l95 = (median_fold_enrichment - quantile(first_fold_enrichments,0.025))/median_fold_enrichment
lower_95 = adjusted_fold_enrichment - (adjusted_fold_enrichment*l95)

u95 = (quantile(first_fold_enrichments,0.975) - median_fold_enrichment)/median_fold_enrichment
upper_95 = adjusted_fold_enrichment + (adjusted_fold_enrichment*u95)

p = length(which(observed_fold_enrichment<permute))/length(permute)
fold = c(adjusted_fold_enrichment,lower_95,upper_95,p)
write.table(fold,argv$out)