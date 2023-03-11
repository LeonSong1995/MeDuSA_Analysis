#! /usr/lib/R/bin/Rscript
library(SnapATAC)
library(viridisLite)
library(ggplot2)
library(parallel)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(patchwork)
library(GenomicRanges)
library(dplyr)
library(ggforce)

snap_file_list = c('../ESO_ATAC.snap')

sample = sapply(snap_file_list,function(snap_file){strsplit(strsplit(snap_file,split='/')[[1]][8],split='.snap')[[1]][1]})
sample = sapply(sample,function(i){paste(strsplit(i,split='_')[[1]][-1],collapse='_')})
names(sample)=NULL
sample="esophagus_mucosa_SM-A9HOR"
x.sp = createSnap(file=snap_file_list,sample=sample,num.cores=4)
sample = paste(sample,1,sep = '_')

#Barcode selection
barcodes = as.data.frame(fread("../GSE184462_metadata.tsv",head=TRUE,sep='\t'))
barcodes = barcodes[2:nrow(barcodes),]
barcodes=barcodes[barcodes$sample %in% sample,]
barcodes$cellID = sapply(barcodes$cellID,function(i){strsplit(i,split='[+]')[[1]][2]})
rownames(barcodes) = barcodes$cellID
x.sp = x.sp[which(x.sp@barcode %in% barcodes$cellID),]
x.sp@metaData = barcodes[x.sp@barcode,]
#Add cell-by-bin matrix
x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=4)
x.sp = makeBinary(x.sp, mat="bmat")

#Bin filtering
black_list = as.data.frame(fread("../hg38-blacklist.v2.bed"))
black_list.gr = GRanges(black_list[,1],IRanges(black_list[,2], black_list[,3]))
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr))
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]}
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))]
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature)
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]}
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1)
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)
idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
x.sp = x.sp[, idy, mat="bmat"]


#Dimensionality reduction
x.sp = runDiffusionMaps(obj=x.sp,input.mat="bmat",num.eigs=20)
plotDimReductPW(obj=x.sp, eigs.dims=1:20,point.size=0.3,
  point.color="grey",point.shape=19,point.alpha=0.6,down.sample=5000,
  pdf.file.name=NULL, pdf.height=7, pdf.width=7)


##Graph-based clustering
x.sp = runKNN(obj=x.sp,eigs.dims=1:10,k=15)
x.sp=runCluster(obj=x.sp,tmp.folder=tempdir(),louvain.lib="R-igraph",seed.use=10)
table(x.sp@cluster)
x.sp@metaData$cluster = x.sp@cluster

#Visualization
x.sp = runViz(obj=x.sp, tmp.folder=tempdir(),dims=2,eigs.dims=1:5,
  method="umap",seed.use=10)

draw = data.frame(x.sp@umap,cluster=x.sp@cluster)
dat = aggregate(draw[,1:2],by=list(draw$cluster),FUN=median)
ggplot(draw,aes(x=umap.1,y=umap.2))+
  xlab('umap-1')+
  ylab('umap-2')+
  geom_point(aes(col=cluster),size=0.5)+
  geom_text(data=dat,aes(x=umap.1,y=umap.2),label=dat$Group.1,size=5,check_overlap = T)+
  theme(legend.position = 'none',
        panel.border = element_blank(), axis.line = element_line(size = 0.8))


genes = as.data.frame((fread('../GENCODE_GRCh38_df.txt')))
genes = genes[,-1]
genes = genes %>% filter(feature=="gene")
genes.gr = GRanges(genes[,'seqname'],IRanges(genes[,'start'], genes[,'end']), name=genes[,'gene_name'])
seqlevelsStyle(genes.gr) = 'UCSC'
marker.genes = c("KRT5","SCPEP1","ITGA6","MKI67","KRT4","ECM1")
genes.sel.gr = genes.gr[which(genes.gr$name %in% marker.genes)]



##cell-type annotation (re-add all bins)
x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=4)
x.sp = makeBinary(x.sp, mat="bmat")
x.sp = createGmatFromMat(obj=x.sp, input.mat="bmat",genes=genes.sel.gr,do.par=TRUE,num.cores=4)
x.sp = scaleCountMatrix(obj=x.sp,mat="gmat",method = "RPM",cov=x.sp@metaData$logUMI)
x.sp = runMagic(obj=x.sp,input.mat="gmat",step.size=3)
#visulaize marker gene
draw = data.frame(scale(x.sp@gmat),cluster=x.sp@cluster)
draw = melt(draw)
colnames(draw) = c('cluster','marker','value')
draw$marker = factor(draw$marker,levels = c("KRT5","SCPEP1","ITGA6","MKI67","KRT4","ECM1"))
p1 = ggplot(draw,aes(x=cluster,y=value))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1))+
  geom_violin(aes(fill=marker))+
  facet_wrap(~marker,scales = 'free')+
  theme(legend.position = 'none',
        axis.title.y = element_blank(),
        strip.background.x = element_blank(),
        strip.placement = 'inside',
        panel.border = element_blank(), axis.line = element_line(size = 0.8))

mk = melt(as.data.frame(scale(x.sp@gmat)))
colnames(mk) = c('marker','value')
draw = data.frame(do.call(rbind, replicate(ncol(x.sp@gmat),data.frame(x.sp@umap,cluster=x.sp@cluster), simplify=FALSE)),mk)
draw$marker = factor(draw$marker,levels = c("KRT5","SCPEP1","ITGA6","MKI67","KRT4","ECM1"))
label = aggregate(x.sp@umap,by=list(x.sp@cluster),FUN=median)
label = do.call(rbind, replicate(ncol(x.sp@gmat),label, simplify=FALSE))
label$marker = rep(c("KRT5","SCPEP1","ITGA6","MKI67","KRT4","ECM1"),each=length(unique(x.sp@cluster)))
colnames(label) = c('cluster','x','y','marker')
label$marker =  factor(label$marker,levels = c("KRT5","SCPEP1","ITGA6","MKI67","KRT4","ECM1"))
label$group = rep('A',nrow(label))
label$group[label$cluster %in% c(8,11,12,6,3,10,9,17,4)]="B"

p2 = ggplot(draw,aes(x=umap.1,y=umap.2))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1))+
  facet_wrap(~marker,scales = 'free')+
  xlab('umap-1')+
  ylab('umap-2')+
  geom_point(aes(col=value),size=0.1)+
  theme(strip.background.x = element_blank(),
        legend.title = element_blank(),
        strip.placement = 'inside',
        panel.border = element_blank(), axis.line = element_line(size = 0.8))+
  geom_text(data=label,mapping=aes(x=x,y=y,label=cluster),hjust=0.5,check_overlap = T)+
  scale_color_gradientn(colours = terrain.colors(20))
patchwork::wrap_plots(p1,p2,ncol=1)

##select the epi cells (re-analysis)
x.sp = x.sp[x.sp@cluster %in% c(4,16,8,18,3,6,10,9,7),]
#Bin filtering
black_list = as.data.frame(fread("../hg38-blacklist.v2.bed"))
black_list.gr = GRanges(black_list[,1],IRanges(black_list[,2], black_list[,3]))
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr))
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]}
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))]
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature)
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]}
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1)
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)
idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
x.sp = x.sp[, idy, mat="bmat"]


#Dimensional reduction
set.seed(10)
x.sp = runDiffusionMaps(obj=x.sp,input.mat="bmat",num.eigs=50)
plotDimReductPW(obj=x.sp, eigs.dims=1:50,point.size=0.3,
                point.color="#a1d99b",point.shape=19,point.alpha=0.6,down.sample=5000,
                pdf.file.name='../eig_ATAC.png', pdf.height=8, pdf.width=8)

##Graph-based clustering
x.sp = runKNN(obj=x.sp,eigs.dims=1:10,k=15)
x.sp=runCluster(obj=x.sp,tmp.folder=tempdir(),louvain.lib="R-igraph",seed.use=10)
table(x.sp@cluster)
x.sp@metaData$cluster = x.sp@cluster

#cell-trajectory annotation (re-add all bins)
genes = as.data.frame((fread('../GENCODE_GRCh38_df.txt')))
genes = genes[,-1]
genes = genes %>% filter(feature=="gene")
genes.gr = GRanges(genes[,'seqname'],IRanges(genes[,'start'], genes[,'end']), name=genes[,'gene_name'])
seqlevelsStyle(genes.gr) = 'UCSC'
marker.genes = c("KRT5","KRT4","ECM1")
genes.sel.gr = genes.gr[which(genes.gr$name %in% marker.genes)]
x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=4)
x.sp = makeBinary(x.sp, mat="bmat")
x.sp = createGmatFromMat(obj=x.sp, input.mat="bmat",genes=genes.sel.gr,do.par=TRUE,num.cores=4)
x.sp = scaleCountMatrix(obj=x.sp,mat="gmat",method = "RPM",cov=x.sp@metaData$logUMI)
x.sp = runMagic(obj=x.sp,input.mat="gmat",step.size=3)

#pseudo-time analyses
library(slingshot)
exp = as.matrix(x.sp@gmat)
exp = melt(exp)[,2:3]
space = data.frame(pc1=x.sp@smat@dmat[,1],pc2=x.sp@smat@dmat[,2])
ggplot(space,aes(x=pc1,y=pc2))+
  geom_point(aes(col=x.sp@cluster))



sds = getLineages(space, x.sp@cluster, start.clus = '2')
sds = getCurves(sds)
path = as.data.frame(sds@metadata$curves$Lineage1$s)
Pseudotime = sds@metadata$curves$Lineage1$lambda
Pseudotime = - Pseudotime
Pseudotime = Pseudotime+abs(min(Pseudotime))
col=c("#ffff99","#edf8b1","#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8")
col = col[seq(length(col),1,-1)]
p1 = ggplot(as.data.frame(space),aes(x=pc1,y=pc2))+
  xlab('PATH1')+
  ylab('PATH2')+
  geom_point(size=2,aes(col=Pseudotime))+
  theme(panel.border = element_blank(), axis.line = element_line(size = 0.8),
        legend.title = element_text(size = 10),
        legend.position = c(0.1,0.15),
        axis.title.y = element_text(size = 12))+
  geom_path(data = path,aes(x=pc1,y=pc2),size=1)+
  annotate('text',x=-0.02,y=0.03,label='Basal layer',size=5)+
  annotate('text',x=0.03,y=0.045,label='Outter layer',size=5)+
  scale_color_gradientn(colours = col)
saveRDS(summary,'../pmat_time.rds')


ct = as.vector(x.sp@cluster)
KRT5 = data.frame('Pseudotime'=Pseudotime,'value'=scale(x.sp@gmat[,'KRT5']),'ct'=ct,'g'=rep('KRT5',length(ct)))
KRT4 = data.frame('Pseudotime'=Pseudotime,'value'=scale(x.sp@gmat[,'KRT4']),'ct'=ct,'g'=rep('KRT4',length(ct)))
ECM1 = data.frame('Pseudotime'=Pseudotime,'value'=scale(x.sp@gmat[,'ECM1']),'ct'=ct,'g'=rep('ECM1',length(ct)))

draw = rbind(KRT5,KRT4,ECM1)
draw$g = factor(draw$g,levels = c("KRT5","KRT4","ECM1"))
p2 = ggplot(draw,aes(x=-Pseudotime,y=value))+
  facet_wrap(~g,ncol=1,scales = 'free')+
  geom_point(aes(col=g),size=0.5)+
  theme(panel.border = element_blank(), axis.line = element_line(size = 0.8),
        strip.background = element_blank(),
        legend.position = 'none',
        strip.placement = "inside",
        legend.title = element_text(size = 10),
        axis.title.y = element_blank())+
  xlab('Pseudotime')+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  scale_color_manual(values =c("#66c2a5","#fc8d62","#8da0cb"))+
  geom_smooth(method='loess',se = T,col='black',linetype='dashed',size=0.8)

patchwork::wrap_plots(p1,p2,ncol=2,widths = c(2,1))

ggsave(filename = "../pseudotime.pdf",
       device = NULL,width = 14,height = 10)


## Call peaks
clusters.sel = unique(x.sp@metaData$cluster)
tempPATH = strsplit(snap_file_list,split='.snap')[[1]][1]
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i])
  runMACS(
    obj=x.sp[which(x.sp@metaData$cluster==clusters.sel[i]),],
    output.prefix=paste0(tempPATH, gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="../snaptools",
    path.to.macs="../macs2",
    gsize="hs",
    buffer.size=500,
    num.cores=2,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir()
  )
}, mc.cores=2)


## add the peak
peaks.names = paste0(paste0(tempPATH,clusters.sel),'_peaks.narrowPeak')
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = reduce(Reduce(c, peak.gr.ls))
peaks.df = as.data.frame(peak.gr)[,1:3]

write.table(peaks.df,file = paste0(tempPATH,".peaks.combined.bed"),append=FALSE,
            quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".",
            row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

path.to.snaptools="../snaptools"
system2(command = path.to.snaptools,
        args = c("snap-add-pmat",
                 "--snap-file", snap_file_list,
                 "--peak-file", paste0(tempPATH,".peaks.combined.bed")))



## peak specificity annotation
x.sp = addPmatToSnap(x.sp, num.cores=2)
index = (colSums(x.sp@pmat)/nrow(x.sp@pmat))>0.01
pmat = x.sp@pmat[,index]
saveRDS(x.sp,'../x.sp_epi.rds')
peak = x.sp@peak
pmat = x.sp@pmat

score = sapply(seq(1:ncol(pmat)),function(i){
  print(i)
  y = pmat[,i]
  y = qnorm((rank(y, na.last="keep")-0.5)/sum(!is.na(y)))
  coef = (summary.gam(mgcv::gam(y~s(Pseudotime,fx = F))))
  coef$chi.sq
})

peak$score = rep(0,length(peak$name))
peak$p = pchisq(score,df = 1,lower.tail = F)
peak$score[index] = score[,3]
peak$p[index] = score[,4]
saveRDS(peak,'/peak.rds')








