#Integrate datasets of six species,taking EXs as example
print(paste("Start time:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(harmony)
library(RColorBrewer)
library(ggpubr)
library(data.table)
library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 100000 * 1024^5)

source("functions/function_merge.R")
source('functions/function_integ.R')

readobj <- function(inrds,inmeta,sub=NULL,incol,invert=F,lab=NULL){
obj <- readRDS(inrds)
meta <- read.table(inmeta,sep='\t',header=T)
obj <- subset(obj,cells=rownames(meta))
obj@meta.data <- meta
if (!is.null(sub)) {
c1 <- which(meta[,incol] %in% sub)
if (invert) {
obj <- subset(obj,cells=rownames(meta)[-c1])
}else{
obj <- subset(obj,cells=rownames(meta)[c1])
}
}
if (!is.null(lab)) {
obj$lab <- lab
}
Idents(obj) <- 'sub'
obj <- subset(obj,downsample = 1000)
return(obj)
}

ob1 <- readobj(inrds='03.human_homology_rds/ma_human_homology.rds',
               inmeta='01.metadata/ma_metadata.xls',
               sub=c('EX'),incol='major',invert=F,lab='ma')

ob2 <- readobj(inrds='03.human_homology_rds/hs.rds',
               inmeta='01.metadata/hs_metadata.xls',
               sub=c('EX'),incol='major',invert=F,lab='hs')

ob3 <- readobj(inrds='03.human_homology_rds/mm_human_homology.rds',
               inmeta='01.metadata/mm_metadata.xls',
               sub=c('EX'),incol='major',invert=F,lab='mm')

ob4 <- readobj(inrds='03.human_homology_rds/zf_human_homology.rds',
               inmeta='01.metadata/zf_metadata.xls',
               sub=c('EX'),incol='major',invert=F,lab='zf')

ob5 <- readobj(inrds='03.human_homology_rds/tt-n_human_homology.rds',
               inmeta='01.metadata/tt_neuron_metadata.xls',
               sub=c('EX'),incol='major',invert=F,lab='tt')

ob6 <- readobj(inrds='03.human_homology_rds/li-n_human_homology.rds',
               inmeta='01.metadata/li_neuron_metadata.xls',
               sub=c('EX'),incol='major',invert=F,lab='li')

gene1 <- Reduce(intersect,list(rownames(ob1),rownames(ob2),rownames(ob3), rownames(ob4),rownames(ob5),rownames(ob6)))
all <- merge(ob1,c(ob2,ob3,ob4,ob5,ob6))
all <- subset(all,features=gene1)
all$library <- all$lab

all <- sct_integ(all,group.by="lab",is.mt=F,is.subset=F)
dassay <- 'RNA'
n <- length(unique(all$library))
tmp <- all
all@assays$integrated@scale.data<-as.matrix(all@assays$integrated@scale.data[1:2,])
saveRDS(all, paste0(n, "libs_integ.rds"))
matx <- all@meta.data
write.table(matx, paste0(n, "sets_cell_info.xls"), sep="\t", quote=F)