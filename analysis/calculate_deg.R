library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 100000 * 1024^5)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(stringr)
library(cowplot)
library(Seurat)
library(dplyr)
source('functions/function_deg.R')
#Calculate DEGs using Seurat
blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)
plot_deg <- function(obj,ident='Celltype',only.pos=T,min.pct = 0.25,logfc.threshold = 0.5,
        assay='RNA',label='test',ntop=5,scale=T,inlevels=NULL,sub='transcript',lname=NULL,invert=F) {
if (!is.null(lname)) {
c1 <- which(obj@meta.data[,ident] %in% lname)
if (invert) {
obj <- subset(obj,cells=Cells(obj)[-c1])
}else{
obj <- subset(obj,cells=Cells(obj)[c1])
}
}
Idents(obj) <- ident
if (is.null(inlevels)) {
inlevels <- sort(levels(obj))
}
obj@meta.data[,ident] <- factor(obj@meta.data[,ident],levels=inlevels)
Idents(obj) <- ident

if (dim(obj)[2] > 50000) {
obj <- subset(obj,downsample = 5000)
}
write.table(obj@meta.data, paste0(label,'_subset_metadata.xls'),sep='\t',quote=F)

DefaultAssay(obj) <- assay
obj <- NormalizeData(obj)
if (scale) {
obj <- ScaleData(obj,features=rownames(obj))
}

all_markers <- FindAllMarkers(obj,only.pos=only.pos,min.pct = min.pct,logfc.threshold = logfc.threshold,verbose = F)

all_markers$FC <- 2^(all_markers$avg_log2FC)
tmp <- all_markers
bklev <- levels(obj)
write.table(all_markers,paste0(label,'_deg_markers.xls'),sep='\t',quote=F)
try({
c1 <- grep(sub,all_markers$gene)
if (length(c1)!=0) {
all_markers <- all_markers[-c1,]
}

all_markers$cluster <- factor(all_markers$cluster,levels=inlevels)

all_markers <- all_markers[order(all_markers$cluster),]
if (!is.null(ntop)) {
top5 <- all_markers %>% group_by(cluster) %>% top_n(n=ntop, wt=avg_log2FC)
}else{
top5 <- all_markers
}
blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)
obj@active.ident <- factor(obj@active.ident,levels=inlevels)

hei <- ceiling(length(unique(top5$gene))*0.2)
pdf(paste0(label,'_deg_dotplot.pdf'),10,hei)
print(
DotPlot(obj,features=unique(top5$gene),assay=assay)+ coord_flip()+
        scale_color_gradientn(colours = blu)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

)
dev.off()
})
all_markers <- tmp

obj@active.ident <- factor(obj@active.ident,levels=bklev)
c1 <- grep(sub,all_markers$gene)
if (length(c1)!=0) {
all_markers <- all_markers[-c1,]
}

if (!is.null(ntop)) {
top5 <- all_markers %>% group_by(cluster) %>% top_n(n=ntop, wt=avg_log2FC)
}else{
top5 <- all_markers
}

hei <- ceiling(length(unique(top5$gene))*0.16)
blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)

jpeg(paste0(label,'_deg_heatmap.jpg'),10*100,hei*100)
print(DoHeatmap(obj,features=unique(top5$gene))+scale_fill_gradientn(colours = blu))
dev.off()
write.table(top5,paste0(label,'_figures_gene_use.xls'),sep='\t',quote=F)

obj$test <- 'test'
Idents(obj) <- 'test'
pdf(paste0(label,'_use_gene_list.pdf'),10,hei)
print(DoHeatmap(obj,features=unique(top5$gene),cells=sample(Cells(obj),100))+scale_fill_gradientn(colours = blu))
dev.off()

}

ntop <- 10
try({
all <- readRDS('06.pfc/macaca_pfc_20230215_diet.rds')
meta <- read.table('01.rds/all_metadata_20230821.xls',sep='\t',header=T)
meta=subset(meta,order=='2macaque')
meta$major=meta$Sub
all=subset(all,cells=rownames(meta))
all@meta.data <- meta
lname <- NULL
invert=T
inlevels <- NULL
plot_deg(all,ident='major',only.pos=T, min.pct = 0.1,logfc.threshold = 0.1, assay='RNA',label='macaca',ntop=ntop,scale=T,inlevels=inlevels,invert=invert,lname=lname)
})

try({
all <- readRDS('06.pfc/human_ref_diet.rds')
meta <- read.table('01.rds/all_metadata_20230821.xls',sep='\t',header=T)
meta=subset(meta,order=='1human')
meta$major=meta$Sub
all=subset(all,cells=rownames(meta))
all@meta.data <- meta
lname <- NULL
invert=T
inlevels <- NULL
plot_deg(all,ident='major',only.pos=T, min.pct = 0.1,logfc.threshold = 0.1, assay='RNA',label='human',ntop=ntop,scale=T,inlevels=inlevels,invert=invert,lname=lname)
})


try({
all <- readRDS('06.pfc/mouse_allen.rds')
meta <- read.table('01.rds/all_metadata_20230821.xls',sep='\t',header=T)
meta=subset(meta,order=='3mouse')
meta$major=meta$Sub
all=subset(all,cells=rownames(meta))
all@meta.data <- meta
lname <- NULL
invert=T
inlevels <- NULL
plot_deg(all,ident='major',only.pos=T, min.pct = 0.1,logfc.threshold = 0.1, assay='RNA',label='mouse',ntop=ntop,scale=T,inlevels=inlevels,invert=invert,lname=lname)

})


if (TRUE) {
try({
all <- readRDS('06.pfc/lizard.rds')
meta <- read.table('01.rds/all_metadata_20230821.xls',sep='\t',header=T)
meta=subset(meta,order=='6lizard')
meta$major=meta$Sub
all=subset(all,cells=rownames(meta))
all@meta.data <- meta
lname <- c('NPC')
invert=T
inlevels <- NULL
plot_deg(all,ident='major',only.pos=T, min.pct = 0.1,logfc.threshold = 0.1, assay='RNA',label='lizard',ntop=ntop,scale=T,inlevels=inlevels,invert=invert,lname=lname)
})
try({
all <- readRDS('06.pfc/turtle.rds')
meta <- read.table('01.rds/all_metadata_20230821.xls',sep='\t',header=T)
meta=subset(meta,order=='5turtle')
meta$major=meta$Sub
all=subset(all,cells=rownames(meta))
all@meta.data <- meta
lname <- c('NPC','doublets')
invert=T
inlevels <- NULL
plot_deg(all,ident='major',only.pos=T, min.pct = 0.1,logfc.threshold = 0.1, assay='RNA',label='turtle',ntop=ntop,scale=T,inlevels=inlevels,invert=invert,lname=lname)
})

try({
all <- readRDS('06.pfc/bird_zf.rds')
meta <- read.table('01.rds/all_metadata_20230821.xls',sep='\t',header=T)
meta=subset(meta,order=='4zebrafinch')
meta$major=meta$Sub
all=subset(all,cells=rownames(meta))
all@meta.data <- meta
lname <- c('NPC')
invert=T
inlevels <- NULL
plot_deg(all,ident='major',only.pos=T, min.pct = 0.1,logfc.threshold = 0.1, assay='RNA',label='zebrafinch',ntop=ntop,scale=T,inlevels=inlevels,invert=invert,lname=lname)

})
}
try({
all <- readRDS('06.pfc/02.analysis/03.amniote_homology/03.cluster/01.EX/6libs_integ.rds')
meta <- read.table('06.pfc/02.analysis/03.amniote_homology/03.cluster/01.EX/6sets_cell_info.xls',sep='\t',header=T,comment='')
all@meta.data <- meta
lname <- NULL
all$lab_major <- paste0(all$lab,'_',all$major)
invert=T
inlevels <- NULL
plot_deg(all,ident='lab_major',only.pos=T, min.pct = 0.1,logfc.threshold = 0.1, assay='RNA',label='ex',ntop=ntop,scale=T,inlevels=inlevels,invert=invert,lname=lname)

})

try({
all <- readRDS('06.pfc/02.analysis/03.amniote_homology/03.cluster/02.IN/6libs_integ.rds')
meta <- read.table('06.pfc/02.analysis/03.amniote_homology/03.cluster/02.IN/6sets_cell_info.xls',sep='\t',header=T,comment='')
all@meta.data <- meta
lname <- NULL
all$lab_major <- paste0(all$lab,'_',all$major)
invert=T
inlevels <- NULL
plot_deg(all,ident='lab_major',only.pos=T, min.pct = 0.1,logfc.threshold = 0.1, assay='RNA',label='in',ntop=ntop,scale=T,inlevels=inlevels,invert=invert,lname=lname)
})

#Calculate DEGs using pesudeobulk methods
plot_deg <- function(obj,ident='Celltype',only.pos=T,min.pct = 0.25,logfc.threshold = 0.5,assay='RNA',label='test',ntop=5,scale=T,inlevels=NULL,sub='transcript',lname=NULL,invert=F) {
if (!is.null(lname)) {
c1 <- which(obj@meta.data[,ident] %in% lname)
if (invert) {
obj <- subset(obj,cells=Cells(obj)[-c1])
}else{
obj <- subset(obj,cells=Cells(obj)[c1])
}
}
Idents(obj) <- ident
if (is.null(inlevels)) {
inlevels <- sort(levels(obj))
}

obj@meta.data[,ident] <- factor(obj@meta.data[,ident],levels=inlevels)
Idents(obj) <- ident

if (dim(obj)[2] > 50000) {
obj <- subset(obj,downsample = 5000)
}
write.table(obj@meta.data, paste0(label,'_subset_metadata.xls'),sep='\t',quote=F)

DefaultAssay(obj) <- assay
obj <- NormalizeData(obj)
if (scale) {
obj <- ScaleData(obj,features=rownames(obj))
}


all_markers <- run_de_multi(sc=obj,replicate=NULL,cell_type=NULL,incol=ident,FC=1.5)

tmp <- all_markers
bklev <- levels(obj)
write.table(all_markers,paste0(label,'_deg_markers.xls'),sep='\t',quote=F)
try({
c1 <- grep(sub,all_markers$gene)
if (length(c1)!=0) {
all_markers <- all_markers[-c1,]
}

all_markers$cluster <- factor(all_markers$cluster,levels=inlevels)

all_markers <- all_markers[order(all_markers$cluster),]
if (!is.null(ntop)) {
top5 <- all_markers %>% group_by(cluster) %>% top_n(n=ntop, wt=avg_log2FC)
}else{
top5 <- all_markers
}
blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)
obj@active.ident <- factor(obj@active.ident,levels=inlevels)

hei <- ceiling(length(unique(top5$gene))*0.2)
pdf(paste0(label,'_deg_dotplot.pdf'),10,hei)
print(
DotPlot(obj,features=unique(top5$gene),assay=assay)+ coord_flip()+
        scale_color_gradientn(colours = blu)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

)
dev.off()
})

all_markers <- tmp

obj@active.ident <- factor(obj@active.ident,levels=bklev)
c1 <- grep(sub,all_markers$gene)
if (length(c1)!=0) {
all_markers <- all_markers[-c1,]
}

if (!is.null(ntop)) {
top5 <- all_markers %>% group_by(cluster) %>% top_n(n=ntop, wt=avg_log2FC)
}else{
top5 <- all_markers
}

hei <- ceiling(length(unique(top5$gene))*0.16)
blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)

jpeg(paste0(label,'_deg_heatmap.jpg'),10*100,hei*100)
print(DoHeatmap(obj,features=unique(top5$gene))+scale_fill_gradientn(colours = blu))
dev.off()
write.table(top5,paste0(label,'_figures_gene_use.xls'),sep='\t',quote=F)

obj$test <- 'test'
Idents(obj) <- 'test'
pdf(paste0(label,'_use_gene_list.pdf'),10,hei)
print(DoHeatmap(obj,features=unique(top5$gene),cells=sample(Cells(obj),100))+scale_fill_gradientn(colours = blu))
dev.off()

}

ntop <- 10

try({
all <- readRDS('06.pfc/mouse_allen.rds')
meta <- read.table('01.rds/all_metadata_20230821.xls',sep='\t',header=T)
meta=subset(meta,order=='3mouse')
meta$major=meta$Sub
all=subset(all,cells=rownames(meta))
all@meta.data <- meta
all =subset(all,subset=Sub!='PERI')
lname <- NULL
invert=T
inlevels <- NULL
plot_deg(all,ident='major',only.pos=T, min.pct = 0.1,logfc.threshold = 0.1, assay='RNA',label='mouse',ntop=ntop,scale=T,inlevels=inlevels,invert=invert,lname=lname)
})
