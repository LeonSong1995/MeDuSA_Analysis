#! /usr/lib/R/bin/Rscript
rm(list=ls())
library(slingshot)
library(MeDuSA)
library(Seurat)
library(reshape2)
library(ggplot2)
library(data.table)
library(ggplot2)
library(gaston)
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
sce = readRDS('../Ref_melanoma.rds')
gene = readRDS('../gene_melanoma.rds')


#Load the bulk data of TCGA-
bulk = readRDS('../bulk_melanoma_TCGA.rds')

# Load the Clinical and TCR data 
clinc = readRDS('../Clinical_melanoma_TCGA.rds')
TCR = readRDS("../TCR_melanoma_TCGA.rds")
TCR = TCR[!is.na(TCR$TRB.clones),]
commonID = intersect(colnames(bulk),rownames(TCR))
bulk = bulk[,commonID]
TCR = TCR[commonID,]
TCRC = TCR$TRB.clones;Q3 = quantile(TCRC,0.33);Q6 = quantile(TCRC,0.66)
L = which(TCRC <= Q3)
M = intersect(which(TCRC <= Q6),which(TCRC > Q3))
diseaseCon = rep('High TCR',length(TCRC))
diseaseCon[L] = 'Low TCR'
diseaseCon[M] = 'Medium TCR'

# Run MeDuSA
MeDuSA_obj = MeDuSA(bulk,sce,select.ct='CD8_Tcells',resolution = 50,
                       markerGene= gene, ncpu=4,start=c(1e-5,1e-2),maxiter=1e+4,
                       smooth=TRUE,smoothMethod='loess',span=0.75,neighbor=5,fractional=F)


abundance = as.matrix(MeDuSA_CAR_NS@Estimation$cell_state_abundance)
gene = MeDuSA_CAR_NS@Estimation$markerGene
bmed = MeDuSA_CAR_NS@Estimation$TimeBin
MANOVA_Pro(MeDuSA_CAR_NS,degree = 2,condition = diseaseCon)


draw  = data.frame('ab' = c(abundance),
                   'bin' = rep(seq(1,length(bmed)),ncol(abundance)),
                   'type' = rep(diseaseCon,each=nrow(abundance)))
draw$type = factor(draw$type,levels = c("High TCR","Medium TCR","Low TCR"))
ggplot(mctd,aes(x=bin,y=ab))+
  geom_errorbar(aes(col=type),stat='summary',size=0.5,width=0.5)+
  geom_line(aes(col=type),stat='summary')+
  scale_color_manual(values =  c("#f46d43","#fee090","#abd9e9"),name='disease state')+
  scale_fill_manual(values = c("#f46d43","#fee090","#abd9e9"),name='disease state')
print(p1)

# Survival analysis 
library("survival")
library("survminer")
clinc_SKCM = clinc
clinc_SKCM = clinc_SKCM[!is.na(clinc_SKCM$OS.time),]
exhauState = colMeans(abundance[ceiling(2*nrow(abundance)/3):nrow(abundance),])
draw = clinc_SKCM[intersect(rownames(clinc_SKCM),names(exhauState)),]
draw$exhau = exhauState[match(rownames(draw),names(exhauState))]
draw$state = rep('low',nrow(draw))
Q5 =  quantile(draw$exhau,0.5)
summary.aov(lm(draw$exhau~draw$tumor_status))
aggregate(draw$exhau,by=list(draw$tumor_status),FUN=mean)
draw$state[draw$exhau > Q5] = 'high'
fit = survfit(Surv(draw$OS.time_imp/365,draw$OS_imp) ~ draw$state,data=draw)
HR = summary(coxph(Surv(draw$OS.time_imp/365,draw$OS_imp) ~ state,data=draw))
J = ggsurvplot(fit,pval = F,conf.int = F,xlab = 'Survival time (year)',ylab = 'Overall survival (%)',
               size = 2,palette = c("#fc8d59", "#91bfdb"),legend =  c(0.7,0.85),
               legend.labs = c("High (>50%)", "Low (<50%)"),
               legend.title = 'Exhaustion State',
               ggtheme = theme(
                 plot.title=element_text(family="Arial", size=22, colour="black"),
                 axis.title.x=element_text(family="Arial", size=20,  colour="black"),
                 axis.title.y=element_text(family="Arial", size=20, colour="black"),
                 axis.text.x=element_text(family="Arial", size=20, colour="black"),
                 axis.text.y=element_text(family="Arial", size=25, colour="black"),
                 axis.ticks=element_line(colour="black"),
                 panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 panel.background=element_blank(),
                 axis.line=element_line(size=1),
                 legend.text=element_text(family="Arial",colour="black",size=20),
                 legend.title=element_text(family="Arial",colour="black",size=20),
                 strip.text.x = element_text(family="Arial",size = 15, colour = "black"),
                 strip.text.y = element_text(family="Arial",size = 15, colour = "black"),
                 strip.background = element_rect(color = "white")))
print(J)
p = survdiff(Surv(draw$OS.time/365,draw$OS) ~ state,data=draw)
p


# ICB analysis 
cohort =as.data.frame(readRDS('../Cohort_melanoma_ICB.rds'))
rownames(cohort) = cohort$V1
cohort = cohort[,-1]
cohort = cohort[-c(which(cohort$BR=='MR'),which(cohort$BR=='SD')),]
cohort = cohort[which(cohort$`biopsy site`=='skin'),]
comm = intersect(rownames(cohort),SampleName)
bulk = bulk[,comm]
cohort = cohort[comm,]

sce$cell_type[sce$cell_type!="CD8_Tcells"]='other'
MeDuSA_CAR_NS = MeDuSA(bulk,sce,select.ct='CD8_Tcells',resolution=100,fixCov=NULL,
                       markerGene=NULL,nbins=10,knots=10,family='gaussian',geneNumber=200,
                       CAR=F,phi=c(0.2,0.4,0.6,0.9),method = 'wilcox',
                       ncpu=6,start=c(1e-5,1e-2),maxiter=1e+4,
                       smooth=T,smoothMethod='loess',span=075,neighbor=5,fractional=F)


all = colMeans(abundance[ceiling(2*nrow(abundance)/3):nrow(abundance),])
aa = all[R]
bb = all[NR]
draw = data.frame(c(aa,bb),c(rep('Responders',length(aa)),rep('Progressors',length(bb))))
colnames(draw) = c('ab','type')
p_ICB = ggplot(draw,aes(x=type,y=ab))+
  scale_color_manual(values = c('#1a9850','#f46d43'))+
  geom_boxplot(outlier.color = NA)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.title.y=element_text(family="Arial", size=15, colour="black"))+
  ylab('Cell abundance')+
  geom_point(aes(col=type),position = position_jitter(0.1),size=2)+
  geom_signif(comparisons = list(c("Responders", "Progressors")),
              textsize = 8,test = 't.test',vjust =2,family = 'Arial',y_position = c(0.1))
print(p_ICB)
