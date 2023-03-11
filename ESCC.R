#! /usr/lib/R/bin/Rscript
rm(list=ls())
library(Seurat)
library(MeDuSA)
library(ggplot2)
myfree<-theme_set(theme_bw())
old_theme <- theme_update(
  plot.title = element_text(size=15, face="bold", colour="black"),
  axis.title.x = element_text(size=12,face="bold", colour="black"),
  axis.title.y = element_text(size=12, face="bold",colour="black"),
  axis.text.x = element_text(size=10, colour="black"),
  axis.text.y = element_text(size=10, colour="black"),
  axis.ticks = element_line(colour="black"),
  panel.grid.major = element_blank(),
  plot.margin = unit(c(0.5,0.8,0.2,0.5), "cm"),
  panel.grid.minor = element_blank(),
  legend.text = element_text(colour="black",size=12),
  legend.title = element_text(colour="black",size=12),
  legend.key = element_blank(),
  legend.background = element_blank(),
  strip.background = element_rect(fill="white"),
  strip.text.x = element_text(face="bold",size = 10, colour = "black"),
  strip.text.y = element_text(face="bold",size = 10, colour = "black")
)

#Load the reference data and the marker gene (to fix the result)
sce = readRDS('../Ref_ESCC.rds')
gene = readRDS('../gene_ESCC.rds')


#Load the bulk data of TCGA------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bulk = readRDS('../bulk_TCGA.rds')
meta = readRDS('../meta_TCGA.rds')

#Run MeDuSA
MeDuSA_obj = MeDuSA(bulk,sce,select.ct='Epi',resolution = 50,
                       markerGene= gene, ncpu=4,start=c(1e-5,1e-2),maxiter=1e+4,
                       smooth=TRUE,smoothMethod='loess',span=0.75,neighbor=5,fractional=F)

abundance = as.matrix(MeDuSA_obj@Estimation$cell_state_abundance)
bmed = MeDuSA_obj@Estimation$TimeBin
MANOVA_Pro(MeDuSA_obj,degree = 2,condition = meta)

draw  = data.frame('ab' = c(abundance),
                   'bin' = rep(bmed,ncol(abundance)),
                   'type' = rep(meta,each = nrow(abundance)))


p1 = ggplot(draw ,aes(x=bin,y=ab))+
  theme(legend.title=element_blank(),
        axis.text.x = element_blank(),
        legend.text=element_text())+
  xlab('Cell state')+
  ylab('Cell state abundance')+
  geom_line(aes(col=type),stat='summary',size=2)+
  geom_errorbar(stat='summary',aes(col=type))+
  scale_color_manual(values = c('#fec44f','#6a51a3'),name='disease state')+
  scale_fill_manual(values = c('#fec44f','#6a51a3'),name='disease state')

print(p1)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Load the bulk data of TCGA-GTEx------------------------------------------------------------------------------------
bulk = readRDS('../bulk_TCGA_GTEx.rds')
meta = readRDS('../meta_TCGA_GTEx.rds')

MeDuSA_obj = MeDuSA(bulk,sce,select.ct='Epi',resolution = 50,
                       markerGene= gene, ncpu=4,start=c(1e-5,1e-2),maxiter = 1e+4,
                       smooth=TRUE,smoothMethod='loess',span=0.75,neighbor = 5,fractional=F)

abundance = as.matrix(MeDuSA_obj@Estimation$cell_state_abundance)
bmed = MeDuSA_obj@Estimation$TimeBin
MANOVA_Pro(MeDuSA_obj,degree = 2,condition = meta)

draw  = data.frame('ab' = c(abundance),
                   'bin' = rep(bmed,ncol(abundance)),
                   'type' = rep(meta,each = nrow(abundance)))


p2 = ggplot(draw ,aes(x=bin,y=ab))+
  theme(legend.title=element_blank(),
        axis.text.x = element_blank(),
        legend.text=element_text())+
  xlab('Cell state')+
  ylab('Cell state abundance')+
  geom_line(aes(col=type),stat='summary',size=2)+
  geom_errorbar(stat='summary',aes(col=type))+
  scale_color_manual(values = c('#fec44f','#6a51a3'),name='disease state')+
  scale_fill_manual(values = c('#fec44f','#6a51a3'),name='disease state')

print(p2)
#------------------------------------------------------------------------------------

#Load the bulk data of GEO------------------------------------------------------------------------------------
bulk = readRDS('../bulk_GEO.rds')
meta = readRDS('../meta_GEO.rds')

MeDuSA_obj = MeDuSA(bulk,sce,select.ct='Epi',resolution = 50,
                       markerGene= gene, ncpu=4,start=c(1e-5,1e-2),maxiter = 1e+4,
                       smooth=TRUE,smoothMethod='loess',span=0.75,neighbor = 5,fractional=F)

abundance = as.matrix(MeDuSA_obj@Estimation$cell_state_abundance)
bmed = MeDuSA_obj@Estimation$TimeBin
MANOVA_Pro(MeDuSA_obj,degree = 2,condition = meta)


draw  = data.frame('ab' = c(abundance),
                   'bin' = rep(bmed,ncol(abundance)),
                   'type' = rep(meta,each = nrow(abundance)))


p3 = ggplot(draw ,aes(x=bin,y=ab))+
  theme(legend.title=element_blank(),
        axis.text.x = element_blank(),
        legend.text=element_text())+
  xlab('Cell state')+
  ylab('Cell state abundance')+
  geom_line(aes(col=type),stat='summary',size=2)+
  geom_errorbar(stat='summary',aes(col=type))+
  scale_color_manual(values = c('#fec44f','#6a51a3'),name='disease state')+
  scale_fill_manual(values = c('#fec44f','#6a51a3'),name='disease state')

print(p3)
#------------------------------------------------------------------------------------
