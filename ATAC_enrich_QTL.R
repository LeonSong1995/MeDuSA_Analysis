#! /usr/lib/R/bin/Rscript
library(viridisLite)
library(ggplot2)
library(parallel)
library(GenomicRanges)
library(data.table)
library(argparser)

p = arg_parser("Prepare SMR esd")
p = add_argument(p, "ieQTL", help="")
p = add_argument(p, "expand", help="")
p = add_argument(p, "m", help="")
p = add_argument(p, "out", help="")
p = add_argument(p, "peak", help="")
argv = parse_args(p)


m = as.numeric(argv$m)
expand = as.numeric(argv$expand)


#null: 1 is null o is eQTL
confounders = as.data.frame(readRDS("../Esophagus_Mucosa.Confounders_Table_small.rds"))
pheno = as.data.frame(read.table('../pheno.bed'))


#peaks
peak = readRDS(argv$peak)

#confounders overlaping with peaks
confounders.g = GRanges(seqnames = paste0(paste('b',confounders$CHR,sep="'"),"'"),
               	   ranges = IRanges(start = confounders$POS-expand,end = confounders$POS+expand))
confounders.use = confounders[na.omit(findOverlaps(peak, confounders.g, select="arbitrary")),]


#prase the ieQTL table
ieQTL = as.data.frame(read.table(argv$ieQTL))
colnames(ieQTL) = c('chr','pos','end','id','gene','strand')
info = confounders[match(ieQTL$id,confounders$SNP),]
ieQTL$AF = info$AF
ieQTL$dis_Tss = info$dis_Tss

#score of ieQTL
csd_eQTL = GRanges(seqnames = paste0(paste('b',ieQTL$chr,sep="'"),"'"),
               	   ranges = IRanges(start = ieQTL$pos-expand,end = ieQTL$pos+expand))
x = peak$score[na.omit(findOverlaps(csd_eQTL, peak, select="arbitrary"))]
x = mean(x)

ieQTL  = ieQTL[na.omit(findOverlaps(peak,csd_eQTL, select="arbitrary")),]


#Run permutations
confounders.use = confounders.use[!confounders.use$SNP %in% ieQTL$id,]
subset = sapply(1:nrow(ieQTL),function(i){
	print(i)
	temp = ieQTL[i,]
	eligible = 	intersect(which(temp$AF==confounders.use$AF),which(temp$dis_Tss==confounders.use$dis_Tss))
	list(eligible)
})
subset = subset[lapply(subset,length)!=0]



y = sapply(1:m,function(i){
	print(i)
	id = c(sapply(subset,function(i){sample(i,1)}))
	pseudo = confounders.use[id,]
	pseudoQTL = GRanges(seqnames = paste0(paste('b',pseudo$CHR,sep="'"),"'"),
	               ranges = IRanges(start = pseudo$POS-expand,end = pseudo$POS+expand))
	score = peak$score[na.omit(findOverlaps(pseudoQTL, peak, select="arbitrary"))]
	mean(score)
})
y = na.omit(y)


y_mean = mean(y)
v_y = var(y)
OR = x/y_mean
a = (x/y_mean)^2
b = v_y / (x^2)
c = v_y / (length(y)*y_mean^2)
SE = sqrt(a*(b+c))
lower = OR-1.96*SE
upper = OR+1.96*SE
Delta.fold = c(lower,OR,upper,x,SE)
names(Delta.fold) = c("lower","OR","upper","percent","SE|p")
Delta.fold
write.table(t(Delta.fold),argv$out,row.names=F)
