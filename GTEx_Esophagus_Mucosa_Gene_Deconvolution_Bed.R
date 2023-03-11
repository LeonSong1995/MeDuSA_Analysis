#! /usr/lib/R/bin/Rscript
library(data.table)
library(MeDuSA)
library(dplyr)
source('../prepare_bed_format_tensorQTL.R')

##STEPI:
##Aim: Load the GTEx RNA-seq data and select samples collected from the Esophagus_Mucosa.
##Note: We only use samples from the genotyped individuals.
#Load the Meta data of samples
SampleMap = as.data.frame(fread('../GTEx/sample_discription'))
SampleMap = SampleMap[grep('Esophagus - Mucosa',SampleMap$SMTSD),]
#Load the TPM matritx
TPM = fread('../GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct')
TPMGene_ID = c(TPM[,1])
TPMGene_Name = c(TPM[,2])
TPM = as.data.frame(TPM[,colnames(TPM)[colnames(TPM) %in% SampleMap$SAMPID],with=F])
#Load the Count matritx
Count = fread('../GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct')
CountGene_ID = c(Count[,1])
CountGene_Name = c(Count[,2])
Count = as.data.frame(Count[,colnames(Count) %in% SampleMap$SAMPID,with=F])
#Ensure the sample order is consistent between TPM and Count
sampleID_RNA = intersect(colnames(TPM),colnames(Count))
TPM = TPM[,sampleID_RNA]
Count = Count[,sampleID_RNA]
#Load names of the genotyped individual
fam = as.data.frame(fread('../GTEx_Analysis_838Indiv.fam'))[,1]
#Extract the sample ID of RNA-seq data
sampleID = sapply(sampleID_RNA,function(i){
  				  temp = strsplit(i,split = '-')[[1]]
         		  paste(temp[1],temp[2],sep = '-')})
#Select the genotypes samples and rename the RNA-seq data with sample ID
commonID = intersect(fam,sampleID)
colnames(TPM) = colnames(Count) = sampleID
TPM = TPM[,commonID]; Count = Count[,commonID] ##497 sample left


##STEPII:
##Aim: Run cell-state deconvolution using MeDuSA
##Note: The reference scRNA-seq data use the symbol gene name.
#Name the bulk data using symbol gene name
symbol = as.vector(unlist(TPMGene_Name))
index = !duplicated(symbol)
bulk = TPM[index,];rownames(bulk) = symbol[index]
#Load the scRNA-seq data reference
sce = readRDS('../Ref_Eso.rds')

#Run MeDuSA (using 4 bins)
MeDuSA_obj = MeDuSA(bulk,sce,select.ct='Epi',resolution = 4,
                       markerGene= gene, ncpu=4,start=c(1e-5,1e-2),maxiter=1e+4,
                       smooth=F,smoothMethod='loess',span=0.75,neighbor=5,fractional=F)
abundance = as.matrix(MeDuSA_obj@Estimation$cell_state_abundance)

#RINT normalize and save the abundance with formart required by tensorQTL
for(i in c('bin1','bin2','bin3','bin4')){
	temp = as.data.frame(abundance[i,])
	temp[,1]=qnorm((rank(temp[,1],na.last="keep")-0.5)/sum(!is.na(temp[,1])))
	path = paste0(paste0('..',i),'.abundance')
	write.table(temp,path,quote=F,col.names=F,row.names=T,sep='\t')
}

#Run MeDuSA (sensitivity analyses: using 10 bins)
MeDuSA_obj = MeDuSA(bulk,sce,select.ct='Epi',resolution = 10,
                       markerGene= gene, ncpu=4,start=c(1e-5,1e-2),maxiter=1e+4,
                       smooth=F,smoothMethod='loess',span=0.75,neighbor=5,fractional=F)
abundance = as.matrix(MeDuSA_obj@Estimation$cell_state_abundance)
#RINT normalize and save the abundance with formart required by tensorQTL
for(i in rownames(abundance_10bin)){
	temp = as.data.frame(abundance_10bin[i,])
	temp[,1]=qnorm((rank(temp[,1],na.last="keep")-0.5)/sum(!is.na(temp[,1])))
	path = paste0(paste0('../',i),'.sensitivity.abundance')
	write.table(temp,path,quote=F,col.names=F,row.names=T,sep='\t')
}


##STEPIII:
##Aim: Perform QC of the RNA-seq data and convert them to the tensor-bed format
ENSG_ID = c(unlist(TPMGene_ID))
ENSG_ID = sapply(ENSG_ID,function(i){strsplit(i,split = '[.]')[[1]][1]})
index = !duplicated(ENSG_ID)
TPM  = TPM[index,]; Count  = Count[index,]; ENSG_ID = ENSG_ID[index]
rownames(TPM) = rownames(Count) = ENSG_ID
RINT_TPM = ConvertBED(TPM,Count,tpm_threshold=0.1,count_threshold=5,sample_frac_threshold=0.2,BedTemplate)
RINT = RINT_TPM$RINT
TPM = RINT_TPM$TPM
#NOTE! only keep the chromosome of 1-22
chr= as.character(seq(1,22))
RINT  = RINT[RINT[,1] %in% chr,]; TPM  = TPM[TPM[,1] %in% chr,]
write.table(RINT,'../GTEx_Esophagus_Mucosa_tensorQTL_RINT.bed',row.names=F,quote=F,sep='\t')
write.table(TPM,'../GTEx_Esophagus_Mucosa_tensorQTL_TPM.bed',row.names=F,quote=F,sep='\t')