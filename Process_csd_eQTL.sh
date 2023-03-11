####################################################################################################
##README:
##This script aims to
#1) Prepare the covariates for the cell-state-dependent (standar) eQTL analyses.
#	The covariates includes: PEER factors + genotype PCs + covariates provided by the GTEx consortium
#2) RUN csd/standard-eQTL using tensor QTL
####################################################################################################

####################################################################################################
##step I:
#Running PEER factors (in shell):
####################################################################################################
exp="../GTEx_Esophagus_Mucosa_tensorQTL_RINT.bed.gz"
out="../GTEx_Esophagus_Mucosa_PEER.mxt"
Rscript /run_Peer.R ${exp} -o ${out}

##Generating genotype PCs (in shell):
command="plink --bfile ../Esophagus_Mucosa \
--geno 0.01 --maf 0.01 --hwe 0.000001 --mind 0.01 \
--pca \
--out ../GTEx_Esophagus_Mucosa_PC"


####################################################################################################
##step II: Merging covariates (in R)
####################################################################################################
#Load the PEER factor and genotype PCs and GTEx provided covariates
library(data.table)
#genotype PCs
GPC = as.data.frame(fread('../GTEx_Esophagus_Mucosa_PC.eigenvec'))
rownames(GPC) = GPC[,1]; GPC=GPC[,-c(1,2)]
colnames(GPC) = paste0('pc',seq(1,ncol(GPC)))
GPC = GPC[,1:4]
#GTEx covariates
GTEx_cov = as.data.frame(fread('../Esophagus_Mucosa.v8.covariates.txt'))
rownames(GTEx_cov) = GTEx_cov[,1]; GTEx_cov=GTEx_cov[,-1]
GTEx_cov = t(GTEx_cov[,rownames(GPC)])
#PEER factors
PEER = as.data.frame(fread('../GTEx_Esophagus_PEER.mxt'))
rownames(PEER) = PEER[,1]; PEER=PEER[,-1]
PEER = PEER[rownames(GTEx_cov),]
#merge the cov
cov = rbind(t(GPC),t(GTEx_cov),t(PEER))
write.table(cov,'../covariates/GTEx_Esophagus_Mucosa_tensorQTL_COV.mtx',quote=F,row.names=T,col.names=T,sep='\t')


####################################################################################################
##step III: eQTL mapping (tensorQTL in shell)
#NOTE: the bed file should be compressed using gzip
#gzip /storage/yangjianLab/songliyang/MeDuSA/cell-state_QTL/GTEx_Esophagus_Mucosa/gene_counts/*bed
####################################################################################################
#type1: cell-state-dependent eQTL:
celltype=(bin1 bin2 bin3 bin4)
for ct in "${celltype[@]}"
do
	plink_prefix_path="../Esophagus_Mucosa"
	expression_bed="../GTEx_Esophagus_Mucosa_tensorQTL_RINT.bed.gz"
	covariates_file="../covariates/GTEx_Esophagus_Mucosa_tensorQTL_COV.mtx"
	interactions_file="../${ct}.abundance"
	prefix="../GTEx_Esophagus_Mucosa_${ct}"

	python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
	--covariates ${covariates_file} \
	--interaction ${interactions_file} \
	--mode cis_nominal
done

#type2: standard eQTL:
plink_prefix_path="../Esophagus_Mucosa"
expression_bed="../GTEx_Esophagus_Mucosa_tensorQTL_RINT.bed.gz"
covariates_file="../GTEx_Esophagus_Mucosa_tensorQTL_COV.mtx"
prefix="../GTEx_Esophagus_Mucosa_standard"

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
--covariates ${covariates_file} \
--mode cis_nominal

#type3: cell-state-dependent eQTL sensitivity analyses (10 bins):
celltype=(bin1 bin2 bin3 bin4 bin5 bin6 bin7 bin8 bin9 bin10)
for ct in "${celltype[@]}"
do
	plink_prefix_path="../Esophagus_Mucosa"
	expression_bed="../GTEx_Esophagus_Mucosa_tensorQTL_RINT.bed.gz"
	covariates_file="../GTEx_Esophagus_Mucosa_tensorQTL_COV.mtx"
	interactions_file="../${ct}.sensitivity.abundance"
	prefix="../GTEx_Esophagus_Mucosa_${ct}.sensitivity"

	python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
	--covariates ${covariates_file} \
	--interaction ${interactions_file} \
	--mode cis_nominal

done


####################################################################################################
##step IV: Convert the .parquet format to .txt format (in shell)
####################################################################################################
celltype=(bin1 bin2 bin3 bin4 standard)
for ct in "${celltype[@]}"
do
	infile="../GTEx_Esophagus_Mucosa_${ct}.cis_qtl_pairs.{TASK_ID}.parquet"
	outfile="../GTEx_Esophagus_Mucosa_${ct}.cis_qtl_pairs.{TASK_ID}.txt"

	python3 ../SaveAsTxt.py \
	--input $infile \
	--output $outfile
done


