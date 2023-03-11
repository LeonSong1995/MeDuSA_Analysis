library(slingshot)
library(pheatmap)
library(Seurat)
library(gaston)
library(reshape2)
library(ggplot2)
library(umap)
library(Rtsne)
library(slingshot)
library(data.table)

# Load the data
sce = readRDS('../blish_covid.seu.rds')
path = ".."
cm.list = paste0(path, list.files(pattern = "*.matrices.rds", path = path))
cm.files = lapply(cm.list, readRDS)
names(cm.files) = sub(path,"",sub("\\_cell.counts.matrices.rds", "", cm.list))

# pre-processing
cm.pp = mapply(EpicPreHS, cm.files, orig.ident = names(cm.files), SIMPLIFY = F)

# Merge
unlisted <- unlist(cm.pp)

# Merge the emat 
unlisted.sort <- unlisted[grep("emat", names(unlisted))]
args = unlisted.sort
output.name <- EpicJoin(args[[1]], args[[2]], by.x = rownames(args[[1]]),by.y = rownames(args[[2]]), all.x = TRUE, all.y = TRUE)
for(i in seq(3,length(args))){
	print(i)
	output.name <- EpicJoin(output.name, args[[i]], by.x = rownames(output.name),by.y = rownames(args[[i]]), all.x = TRUE, all.y = TRUE)
}
covid_combined.emat <- output.name
colnames(covid_combined.emat) = sapply(colnames(covid_combined.emat),function(i){
	temp = strsplit(i,split='_')[[1]]
	paste(temp[-c(1,length(temp))],collapse='_')
})
colnames(covid_combined.emat)[-grep('HIP',colnames(covid_combined.emat))] = 
paste0('covid_',colnames(covid_combined.emat)[-grep('HIP',colnames(covid_combined.emat))])

saveRDS(covid_combined.emat,'../covid_combined.emat.rds')


# Merge the nmat 
unlisted.sort <- unlisted[grep("nmat", names(unlisted))]
args = unlisted.sort
output.name <- EpicJoin(args[[1]], args[[2]], by.x = rownames(args[[1]]),by.y = rownames(args[[2]]), all.x = TRUE, all.y = TRUE)
for(i in seq(3,length(args))){
	print(i)
	output.name <- EpicJoin(output.name, args[[i]], by.x = rownames(output.name),by.y = rownames(args[[i]]), all.x = TRUE, all.y = TRUE)
}
covid_combined.nmat <- output.name
colnames(covid_combined.nmat) = sapply(colnames(covid_combined.nmat),function(i){
	temp = strsplit(i,split='_')[[1]]
	paste(temp[-c(1,length(temp))],collapse='_')
})
colnames(covid_combined.nmat)[-grep('HIP',colnames(covid_combined.nmat))] = 
paste0('covid_',colnames(covid_combined.nmat)[-grep('HIP',colnames(covid_combined.nmat))])
table(colnames(sce) %in% colnames(covid_combined.nmat))
saveRDS(covid_combined.nmat,'../covid_combined.nmat.rds')

# Merge the smat 
unlisted.sort <- unlisted[grep("smat", names(unlisted))]
args = unlisted.sort
output.name <- EpicJoin(args[[1]], args[[2]], by.x = rownames(args[[1]]),by.y = rownames(args[[2]]), all.x = TRUE, all.y = TRUE)
for(i in seq(3,length(args))){
	print(i)
	output.name <- EpicJoin(output.name, args[[i]], by.x = rownames(output.name),by.y = rownames(args[[i]]), all.x = TRUE, all.y = TRUE)
}
covid_combined.smat <- output.name
colnames(covid_combined.smat) = sapply(colnames(covid_combined.smat),function(i){
	temp = strsplit(i,split='_')[[1]]
	paste(temp[-c(1,length(temp))],collapse='_')
})
colnames(covid_combined.smat)[-grep('HIP',colnames(covid_combined.smat))] = 
paste0('covid_',colnames(covid_combined.smat)[-grep('HIP',colnames(covid_combined.smat))])
table(colnames(sce) %in% colnames(covid_combined.smat))
saveRDS(covid_combined.smat,'../covid_combined.smat.rds')

# Make Seurat object
covid_combined.emat = covid_combined.emat[,colnames(sce)]
covid_combined.nmat = covid_combined.nmat[,colnames(sce)]
covid_combined.smat = covid_combined.smat[,colnames(sce)]
covid_combined <- CreateSeuratObject(counts = covid_combined.emat, min.cells = 10, names.field = 1, names.delim = "\\.")
covid_combined[["spliced"]] <- CreateAssayObject(covid_combined.emat )
covid_combined[["unspliced"]] <- CreateAssayObject(covid_combined.nmat)
covid_combined[["spanning"]] <- CreateAssayObject(covid_combined.smat)
covid_combined$cell.type.fine = sce$cell.type.fine
covid_combined$cell.type.coarse = sce$cell.type.coarse
saveRDS(covid_combined,'../covid_combined.emat.seurat.rds')


Tcell.index = c(grep('CD8 T',covid_combined$cell.type.coarse))
# Reduction
covid_combined_Tcell = covid_combined[,Tcell.index]
covid_combined_Tcell = PercentageFeatureSet(covid_combined_Tcell, pattern = "^MT-", col.name = "percent.mt")
covid_combined_Tcell = PercentageFeatureSet(covid_combined_Tcell, pattern = "^RPS", col.name = "percent.rps")
covid_combined_Tcell = PercentageFeatureSet(covid_combined_Tcell, pattern = "^RPL", col.name = "percent.rpl")
covid_combined_Tcell = PercentageFeatureSet(covid_combined_Tcell, pattern = "^RNA\\d8S5", col.name = "percent.rrna")
covid_combined_Tcell = SCTransform(covid_combined_Tcell, vars.to.regress = c("percent.mt", "percent.rps", "percent.rpl", "percent.rrna", "nCount_RNA", "nFeature_RNA"), verbose = TRUE, return.only.var.genes = TRUE) #expect "iteration limit reached" warning unless suppressed per https://github.com/satijalab/seurat/issues/1426
covid_combined_Tcell = RunPCA(covid_combined_Tcell, verbose = FALSE)
covid_combined_Tcell = RunUMAP(covid_combined_Tcell, dims = 1:20, verbose = FALSE)

# RunVelocity
covid_combined_Tcell = RunVelocity(object = covid_combined_Tcell, deltaT = 1, kCells = 25, fit.quantile = 0.02)
saveRDS(covid_combined_Tcell,'../covid_combined_Tcell.seurat.rds')


covid_combined_Tcell = readRDS('../covid_combined_Tcell.seurat.rds')
ident.colors = (scales::hue_pal())(n = length(x = levels(x = covid_combined_Tcell)))
names(x = ident.colors) = levels(x = covid_combined_Tcell)
cell.colors = ident.colors[Idents(object = covid_combined_Tcell)]
names(x = cell.colors) = colnames(x = covid_combined_Tcell)

umap = Embeddings(object = covid_combined_Tcell, reduction = "umap")
vel = Tool(object = covid_combined_Tcell, slot = "RunVelocity")

pdf(file = "../velocyto.pdf",width = 10,height = 10)
Fig = velocyto.R::show.velocity.on.embedding.cor(emb =umap , vel =vel , n.cores=1,
	n = 200, scale = "sqrt", cell.colors = cell.colors, return.details = T,
    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, 
    min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
    do.par = FALSE, cell.border.alpha = 0.1)
dev.off()


# CD8-cell re-analysis 
## marker genes
TN=c('CCR7','LEF1','TCF7')
TEM=c('GZMK','GZMH', 'GZMM','CXCR4', "CD74", "CXCR5", "CCR4", "CD44")
TEX = c('GZMA', 'GNLY', 'GZMB', 'IFNG','TOX2','HAVCR2')
TEX2 = c('TOX2', 'RUNX2', 'ENTPD1','PDCD1')
TRM=c('ZNF683','ITGAE')
TEMRA=c("GZMH", "GZMB", "GZMM", "GZMA")
TMAIT = c('SLC4A10', 'PRSS35', 'CCR6')
Ty=c('TRDC', 'TRGC2')

CD8 = readRDS('../covid_combined_Tcell.seurat.rds')
out = c('Sex','percent.mt','percent.rps','percent.rpl','percent.rrna','Donor')
CD8 =  SCTransform(CD8, variable.features.n = 2000,vars.to.regress = out,return.only.var.genes = T)
CD8= RunPCA(CD8, features = VariableFeatures(object = CD8))
CD8 = RunUMAP(CD8, dims = 1:20)
CD8 = FindNeighbors(CD8, dims = 1:20)
CD8 = FindClusters(CD8, resolution =0.3)
DimPlot(CD8, reduction = "umap",label = T)
VlnPlot(CD8,features = TN)

#Slingshot 
space = Embeddings(CD8,'umap')
space = as.data.frame(space)
space$ct = as.vector(CD8$SCT_snn_res.0.3[rownames(space)])
sds = getLineages(space[,1:2], space$ct, start.clus = '1')
sds = getCurves(sds)
path = as.data.frame(sds@metadata$curves$Lineage1$s)
time = sds@metadata$curves$Lineage1$lambda
time = time/max(time)

space$ct = as.vector(CD8$SCT_snn_res.0.3[rownames(space)])
ct  = space$ct
ct[ct  %in% c('1')]='CD8_Naive'
ct[ct  %in% c('0','4')]='CD8_EM'
ct[ct  %in% c('3')]='CD8_EM-EX'
ct[ct  %in% c('2')]='CD8_EX'
space$ct =ct
ident_colours = c('#1a9850','#d73027','#2166ac','#fd8d3c')
p1 = ggplot(space,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(col=ct),alpha=0.5)+
  theme(legend.title = element_blank(),
        panel.border = element_blank(), axis.line = element_line(size = 0.8),
        axis.title.y = element_text(size = 12))+
  scale_color_manual(values = ident_colours)+
  geom_path(data =path ,aes(x=UMAP_1,y=UMAP_2),size=1,
            arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")))
print(p1)

TEM = c('HLA-DRB1','CXCR2','CXCR1','CD74','LAG3','IL7R','GZMK','GZMH', 'GZMA','GZMB','CXCR4', "CD74", "CXCR5", "CCR4", "CD44","CTLA4","PDCD1",'CCR7','PCNA','CXCR6')
CD8 = ScaleData(CD8,assay = 'RNA',features = TEM)
plot(time,CD8@assays$RNA@scale.data['CCR7',])
plot(time,CD8@assays$RNA@scale.data['IL7R',])

plot(time,CD8@assays$RNA@scale.data['GZMH',])
plot(time,CD8@assays$RNA@scale.data['HLA-DRB1',])

plot(time,CD8@assays$RNA@scale.data['GZMB',])
plot(time,CD8@assays$RNA@scale.data['GZMK',])

CCR7 = data.frame('pesudos-state'=time,'value'=scale(CD8@assays$RNA@scale.data['CCR7',]),'g'=rep('CCR7',ncol(CD8)))
IL7R = data.frame('pesudos-state'=time,'value'=scale(CD8@assays$RNA@scale.data['IL7R',]),'g'=rep('IL7R',ncol(CD8)))
GZMH = data.frame('pesudos-state'=time,'value'=scale(CD8@assays$RNA@scale.data['GZMH',]),'g'=rep('GZMH',ncol(CD8)))
HLA = data.frame('pesudos-state'=time,'value'=scale(CD8@assays$RNA@scale.data['HLA-DRB1',]),'g'=rep('HLA-DRB1',ncol(CD8)))
GZMB = data.frame('pesudos-state'=time,'value'=scale(CD8@assays$RNA@scale.data['GZMB',]),'g'=rep('GZMB',ncol(CD8)))
PCNA = data.frame('pesudos-state'=time,'value'=scale(CD8@assays$RNA@scale.data['PCNA',]),'g'=rep('PCNA',ncol(CD8)))

draw = rbind(CCR7,IL7R,GZMH,HLA,GZMB,PCNA)
draw$g  = factor(draw$g,levels = c('CCR7','IL7R','GZMH','HLA-DRB1','GZMB','PCNA'))
p2 = ggplot(draw,aes(x=pesudos.state,y=value))+
  facet_wrap(~g,ncol=2,scales = 'free')+
  theme(legend.position = 'none',
        legend.title = element_blank(),
        strip.background = element_blank())+
  xlab('Pseudotime')+
  geom_smooth(method='loess',se = T,size=0.8,aes(col=g,fill=g))+
  ylab('Expression level (Z score)')
print(p2)


CD8$cell_trajectory = time[colnames(CD8)]
CD8$cell_type = 'CD8 T'
saveRDS(CD8,'../scRNA_covid_CD8.rds')

bulk =  readRDS('../covid19_data(1).rds')
marker_gene = MeDuSA::MeDuSA_marker(sce,bulk,method='wilcox',ncpu=4)
saveRDS(CD8,'../gene_COVID.rds')

