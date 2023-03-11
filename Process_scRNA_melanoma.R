#! /usr/lib/R/bin/Rscript
library(slingshot)
library(MeDuSA)
library(Seurat)
library(reshape2)
library(ggplot2)
library(data.table)
library(ggplot2)
library(gaston)


ExhauALL = readRDS('../exhasution_melanoma.rds')

## Load the raw data
raw = read.table('../GSE72056_melanoma_single_cell_revised_v2.txt',header=T)
raw_pro = raw
raw_pro = raw_pro[!duplicated(raw_pro[,1]),]
rownames(raw_pro) = raw_pro[,1]; raw_pro = raw_pro[,-1]
cell_mark = raw_pro[c(1,2,3),]
sample_id  = raw_pro[1,]
raw_pro_eli = raw_pro[-c(1,2,3),]
tpm = (2^raw_pro_eli)-1
tpm = 1e+6*sweep(tpm,2,colSums(tpm),'/')
Tcell = tpm[,colnames(cell_mark[,which(cell_mark[3,]==1)])]
cd4 = colnames(Tcell[,which(Tcell['CD4',]>0)])
cd8 = unique(c(colnames(Tcell[,which(Tcell['CD8A',]>0)]),colnames(Tcell[,which(Tcell['CD8B',]>0)])))
CD4_T  = tpm[,setdiff(cd4,cd8)]
CD8_T  = tpm[,setdiff(cd8,cd4)]
Bcell  =  tpm[,colnames(cell_mark[,which(cell_mark[3,]==2)])]
Macrophages  =  tpm[,names(cell_mark[,which(cell_mark[3,]==3)])]
NK  =  tpm[,tpm(cell_mark[,which(cell_mark[3,]==6)])]
Endothelial  =  tpm[,names(cell_mark[,which(cell_mark[3,]==4)])]
CAF  =  tpm[,names(cell_mark[,which(cell_mark[3,]==5)])]
sce = CreateSeuratObject(tpm)
sce$ct = 'cell'
sce$ct[colnames(sce) %in % colnames(CD8_T)] = 'CD8_Tcells'
saveRDS(sce,'../melanoma_sc.rds')

# Annotate CD8+ T-cell
sce = readRDS('../melanoma_sc.rds')
sce = RunUMAP(sce,dims = 1:20)
CD8 = sce[,Idents(sce)=='CD8_Tcells']
CD8  = NormalizeData(CD8)
CD8 =  FindVariableFeatures(CD8, selection.method = "vst", nfeatures = 500)
CD8 = ScaleData(CD8)
CD8 = RunPCA(CD8, features = VariableFeatures(object = CD8))
CD8 = FindNeighbors(CD8, dims = 1:10)
CD8 = FindClusters(CD8, resolution = 0.2)
CD8 = RunUMAP(CD8, dims = 1:10)
Idents(CD8) = CD8$seurat_clusters
DimPlot(CD8, reduction = "umap")

# Expression Pattern
Naive = c('CCR7','TCF7','LEF1','SELL')
Exhau=ExhauALL[[1]];Exhau2=ExhauALL[[2]];Exhau3=ExhauALL[[3]]
time = colMeans(as.matrix(CD8@assays$RNA@counts[Exhau,])) - colMeans(as.matrix(CD8@assays$RNA@counts[Naive,]))
time2 = colMeans(as.matrix(CD8@assays$RNA@counts[Exhau2,])) - colMeans(as.matrix(CD8@assays$RNA@counts[Naive,]))
time3 = colMeans(as.matrix(CD8@assays$RNA@counts[Exhau3,])) - colMeans(as.matrix(CD8@assays$RNA@counts[Naive,]))

time = time+abs(min(time))
time = time/max(time)
time2 = time2+abs(min(time2))
time2 = time2/max(time2)
time3 = time3+abs(min(time3))
time3 = time3/max(time3)

# Visulaize 
CD8 = ScaleData(CD8,features = c("CCR7","TCF7","SELL","CTLA4","PDCD1","TOX"))
draw.CTL4 = data.frame(time,exp=CD8@assays$RNA@scale.data['CTLA4',],gene=rep('CTLA4',length(time)))
draw.PDCD1 = data.frame(time,exp=CD8@assays$RNA@scale.data['PDCD1',],gene=rep('PDCD1',length(time)))
draw.TOX = data.frame(time,exp=CD8@assays$RNA@scale.data['TOX',],gene=rep('TOX',length(time)))
draw.CCR7 = data.frame(time,exp=CD8@assays$RNA@scale.data['CCR7',],gene=rep('CCR7',length(time)))
draw.TCF7 = data.frame(time,exp=CD8@assays$RNA@scale.data['TCF7',],gene=rep('TCF7',length(time)))
draw.SELL = data.frame(time,exp=CD8@assays$RNA@scale.data['SELL',],gene=rep('SELL',length(time)))
draw = rbind(draw.CTL4,draw.PDCD1,draw.TOX,draw.CCR7,draw.TCF7,draw.SELL)

p1 = ggplot(draw,aes(x=time,y=exp))+
  geom_point(aes(col=gene))+
  geom_smooth()+
  facet_wrap(~gene,scales = 'free_y')

print(p1)  
sce$cell_trajectory = 0
sce$cell_trajectory[colnames(CD8)] = time

bulk = readRDS('../TCGA_melanoma.rds')

# Try the wilcox test, the result is similar. 
gene = MeDuSA::MeDuSA_marker(method = 'gam',sce = CD8,k = 10,family='gaussian',
                     bulk = bulk,geneNumber = 200,ncpu = 4,nbins = 10)
saveRDS(gene,'../gene_melanoma.rds')

