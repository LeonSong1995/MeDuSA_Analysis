#! /usr/lib/R/bin/Rscript
oes = readRDS('../oesophagus_ts.rds')
ct = as.vector(oes$Celltypes_GenomeBiol_2019)
ct[which(ct=='Epi_dividing')]='Epi_basal'
ct[which(ct=='Epi_suprabasal')]='Epi_basal'
Idents(oes) = ct

###bulk(oes)#############
meta = readRDS('../bulk_Sc_meta.rds')
meta = meta[meta[,2]=='Oesophagus',]
meta=meta[meta[,1]%in% unique(oes$donor_time),]
data = readRDS('../tpm.rds')
data = data[,rownames(meta)]
colnames(data) = meta[,1]
data = data[,grep('0h',colnames(data))]


#######ref##############
oes = oes[,oes$donor_time %in% colnames(data)]
oes = NormalizeData(oes, normalization.method = "LogNormalize", scale.factor = 10000)
oes = FindVariableFeatures(oes, selection.method = "vst", nfeatures = 2000)
oes = ScaleData(oes)
oes = RunPCA(oes, features = VariableFeatures(object = oes))
oes = RunUMAP(oes, dims = 1:10)
ump = as.data.frame(Embeddings(oes[["umap"]])[,c(1,2)])
ct = as.vector(Idents(oes))
ct[grep('Epi',ct)] ='Epithelial cell'
ct[grep('T_',ct)] ='T cell'
ct[grep('B_',ct)] ='B cell'
ct[grep('Glands',ct)] ='Glands'
ct[grep('Mono',ct)] ='Monocyte'
ct[grep('Glands',ct)] ='Glands'
ct[grep('vessel',ct)] ='Vessel'
ct[grep('Dendritic_Cells',ct)] ='Dendritic cell'
ct[grep('Mast_cell',ct)] ='Mast cell'
oes$ct = ct

######umap plot for all cells#######################################################################################################################################
col = c('#bebada','#8dd3c7','#fdb462', '#ffffb3','#fb8072', '#80b1d3','#b3de69','#fccde5','#d9d9d9')
ggplot(ump,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(col=ct),size=0.5)+
  theme(legend.position = c(0.1,0.3),
        legend.title = element_blank(),
        panel.border = element_blank(), axis.line = element_line(size = 0.8),
        axis.title.y = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_color_manual(values = col)

VlnPlot(oes,features = c("KRT5","SCPEP1","ITGA6","MKI67","KRT4","ECM1"),group.by = "ct",
        cols = col,pt.size = 0.05,ncol = 3)+
        theme(legend.position = 'none',
              legend.title = element_blank(),
              panel.border = element_blank(), axis.line = element_line(size = 0.8))


######umap plot for epi cells (markers)#######################################################################################################################################
index = grep('Epi',ct)
ggplot(ump[index,],aes(x=UMAP_1,y=UMAP_2))+
  theme(legend.title=element_text(family="Arial",colour="black",size=15))+
  scale_color_gradient(low='#e0e0e0',high='#b2182b',name="KRT5")+
  geom_point(aes(col=oes@assays$RNA@counts['KRT5',index]),alpha=1)

ggplot(ump[index,],aes(x=UMAP_1,y=UMAP_2))+
  theme(legend.title=element_text(family="Arial",colour="black",size=15))+
  scale_color_gradient(low='#e0e0e0',high='#b2182b',name="KRT4")+
  geom_point(aes(col=oes@assays$RNA@counts['KRT4',index]),alpha=1)

ggplot(ump[index,],aes(x=UMAP_1,y=UMAP_2))+
  theme(legend.title=element_text(family="Arial",colour="black",size=15))+
  scale_color_gradient(low='#e0e0e0',high='#b2182b',name="ECM1")+
  geom_point(aes(col=oes@assays$RNA@counts['ECM1',index]),alpha=1)

ct = as.vector(oes$Celltypes_GenomeBiol_2019)[index]
ct[grep('Epi_basal',ct)]='Epithelial basal'
ct[grep('Epi_dividing',ct)]='Epithelial basal'
ct[grep('Epi_stratified',ct)]='Epithelial stratified'
ct[grep('Epi_suprabasal',ct)]='Epithelial basal'
ct[grep('Epi_upper',ct)]='Epithelial upper'

ggplot(ump[index,],aes(x=UMAP_1,y=UMAP_2))+
  guides(color = guide_legend(override.aes = list(size = 10)))+
  theme(legend.title=element_blank(),
        legend.text=element_text(family="Arial",colour="black",size=25))+
  geom_point(aes(col=ct),alpha=1)

#############annotated ct#############
ct = as.vector(oes$Celltypes_GenomeBiol_2019)
ct[grep('Epi_basal',ct)]='Epithelial basal'
ct[grep('Epi_dividing',ct)]='Epithelial basal'
ct[grep('Epi_stratified',ct)]='Epithelial stratified'
ct[grep('Epi_suprabasal',ct)]='Epithelial basal'
ct[grep('Epi_upper',ct)]='Epithelial upper'
oes$ct = ct

##########epi reanalysis#############
oes = readRDS('../oesophagus_part.rds')
sce = oes[,grep('Epi',Idents(oes))]
sce = NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce = FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce = ScaleData(sce)
sce = RunPCA(sce, features = VariableFeatures(object = sce))
sce = RunUMAP(sce,dims = 1:10)
space = Embeddings(sce[["pca"]])
umap = Embeddings(sce[["umap"]])
space[,1] = -space[,1]
df = data.frame('pc1'=space[,1],'pc2'=space[,2],'type'=sce$Celltypes_GenomeBiol_2019)
ct = as.vector(sce$Celltypes_GenomeBiol_2019)
ct[grep('Epi_basal',ct)]='Epithelial basal'
ct[grep('Epi_dividing',ct)]='Epithelial basal'
ct[grep('Epi_stratified',ct)]='Epithelial stratified'
ct[grep('Epi_suprabasal',ct)]='Epithelial basal'
ct[grep('Epi_upper',ct)]='Epithelial upper'


##########cell trajectory#############
sds = getLineages(space, ct, start.clus = 'Epithelial basal')
sds = getCurves(sds)
path = as.data.frame(sds@metadata$curves$Lineage1$s)
Pseudotime =sds@metadata$curves$Lineage1$lambda
Pseudotime = Pseudotime/max(Pseudotime)
col=c("#ffff99","#edf8b1","#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8")
col = col[seq(length(col),1,-1)]
space = as.data.frame(space)
space$Pseudotime = Pseudotime
ggplot(space,aes(x=PC_1,y=PC_2))+
  xlab('PATH-1')+
  ylab('PATH-2')+
  geom_point(aes(col=Pseudotime),size=2)+
  annotate('text',x=-20,y=25,label='Basal layer',size=5)+
  annotate('text',x=20,y=20,label='Outter layer',size=5)+
  theme(strip.background = element_blank(),
        legend.position = c(0.1,0.2),
        panel.border = element_blank(), axis.line = element_line(size = 0.8))+
  geom_path(data = path,aes(x=PC_1,y=PC_2),size=1.5,
            arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")))+
  scale_color_gradientn(colours = col)


##########expression pattern of marker genes#############
KRT5 = data.frame('pesudos-state'=df[,1],'value'=scale(sce@assays$RNA@scale.data['KRT5',]),'g'=rep('KRT5',length(ct)))
KRT4 = data.frame('pesudos-state'=df[,1],'value'=scale(sce@assays$RNA@scale.data['KRT4',]),'g'=rep('KRT4',length(ct)))
ECM1 = data.frame('pesudos-state'=df[,1],'value'=scale(sce@assays$RNA@scale.data['ECM1',]),'g'=rep('ECM1',length(ct)))
draw = rbind(KRT5,KRT4,ECM1)
draw$g  = factor(draw$g,levels = c('KRT5','KRT4','ECM1'))
ggplot(draw,aes(x=pesudos.state,y=value))+
  facet_wrap(~g,ncol=3,scales = 'free')+
  scale_color_gradientn(colours = col)+
  theme(legend.position = 'none',
        legend.title = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(), axis.line = element_line(size = 0.8))+
  geom_point(aes(col=pesudos.state),alpha=1)+
  xlab('Pseudotime')+
  geom_smooth(method='loess',se = F,col='black',linetype='dashed',size=0.8)+
  ylab('Expression level')

