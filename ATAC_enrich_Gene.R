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
peak = readRDS(as.vector(argv$peak))

#null: 1 is null o is eQTL
genes = as.data.frame((fread('../GENCODE_GRCh38_df.txt')))
genes = genes[,-1]
genes = genes[genes$feature=="gene",]

ieQTL = as.data.frame(fread(as.vector(argv$ieQTL)))
if(type=='eQTL'){
	eGene = ieQTL$gene_id[ieQTL$pval_beta<argv$p]
	eGene = sapply(eGene,function(i){strsplit(i,split='[.]')[[1]][1]})
}else if(type=='ieQTL'){
	eGene = ieQTL$gene_id[ieQTL$pval_adj_bh<argv$p]
}

chr = genes$seqname[match(eGene,genes$gene_id)]
start = genes$start[match(eGene,genes$gene_id)]
end = genes$end[match(eGene,genes$gene_id)]



##Run permutations
genesNULL = genes[-match(eGene,genes$gene_id),]
permute = sapply(1:10000,function(i){
	print(i)
	genesRandom = genesNULL[sample(seq(1,nrow(genesNULL)),length(eGene)),]
	pseudoQTL = GRanges(seqnames = paste0(paste('b',genesRandom$seqname,sep="'"),"'"),
	               ranges = IRanges(start = genesRandom$start,end = genesRandom$end))
	r=peak[GenomicRanges::findOverlaps(pseudoQTL,peak)@to,]$score
	mean(r)
})
permute = na.omit(permute)


csd_eQTL = GRanges(seqnames = paste0(paste('b',chr,sep="'"),"'"),
               ranges = IRanges(start = start,end = end))


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