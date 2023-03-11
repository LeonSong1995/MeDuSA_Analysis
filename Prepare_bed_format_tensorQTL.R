library(dplyr)
library(data.table)
message('loading the gtf annotation')
gtf = Rgb::read.gtf('../refdata-gex-GRCh38-2020-A/genes/genes.gtf')
gtf = gtf %>% filter(feature=="gene")
BedTemplate = gtf %>% select('#chr'=seqname,start,end,gene_id)
BedTemplate$gene_id = sapply(BedTemplate$gene_id,function(i){strsplit(i,split = '[.]')[[1]][1]})
BedTemplate$end = BedTemplate$start+1
BedTemplate[,1] = sapply(BedTemplate[,1],function(i){strsplit(i,split="chr")[[1]][2]})
BedTemplate = as.data.frame(BedTemplate)

ConvertBED = function(TPM,Count,tpm_threshold=0.1,count_threshold=6,sample_frac_threshold=0.2,BedTemplate=BedTemplate){
  # expression thresholds
  common = intersect(BedTemplate$gene_id, rownames(TPM))
  TPM = TPM[common,]
  Count = Count[common,]
  message("QC")
  sampleNum = ncol(TPM)
  indexTPM = apply(TPM,1,function(i){length(which(i>=tpm_threshold)) >= sampleNum*sample_frac_threshold})
  indexCount = apply(Count,1,function(i){length(which(i>=count_threshold)) >= sampleNum*sample_frac_threshold})
  eligible = indexTPM & indexCount
  TPM = TPM[eligible,]
  Count = Count[eligible,]

  # apply normalization (TMM and RINT)
  message("TMM and RINT")
  Count = edgeR::DGEList(Count)
  Count = edgeR::calcNormFactors(Count, method = "TMM")
  tmm = edgeR::cpm(Count)
  TMM_RINT = t(apply(tmm,1,function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}))
  TPM_RINT = t(apply(TPM,1,function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}))


  ##write the bed file
  message("Converting to bed")
  Info = BedTemplate[na.omit(match(rownames(TPM),BedTemplate$gene_id)),]
  TPM_RINT = TPM_RINT[Info$gene_id,]
  TMM_RINT = TMM_RINT[Info$gene_id,]
  TPM_RINT = cbind(Info,TPM_RINT)
  TMM_RINT = cbind(Info,TMM_RINT)

  print(dim(TPM_RINT));print(dim(TMM_RINT))
  return(list(TPM_RINT=TPM_RINT,TMM_RINT=TMM_RINT))
}