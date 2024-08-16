
# **********************************************************
# * Author        : HuangFubaoqian
# * Email         : huangbaoqian@genomics.cn
# * Create time   : 2022-05-20 12:37
# * Filename      : function_integ.R
# * Description   : 
# **********************************************************

library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)

library(harmony)
library(cowplot)
library(RColorBrewer)
library(ggpubr)
library(ggsci)
library(data.table)

colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))
len <- 100
cor <- c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"),colpalettes,rainbow(len))
blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)

seurat_integ <- function(obj,group.by=NULL,
                         mt.pattern="^MT-",mt.list=NULL,dim.use=30,mt.cutoff=5,
                         nf.low=500,nf.high=6000,nfeatures=3000,
                         res=1.5,is.mt=NULL,is.subset=T) {
all <- obj
if (is.null(is.mt)) {
if (is.null(mt.list)) {
all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = mt.pattern)
}else{
mt.list <- mt.list[which(mt.list %in% rownames(all))]
all[["percent.mt"]] <- PercentageFeatureSet(all, features = mt.list)
}
}

if (is.subset) {
all <- subset(all, subset = nFeature_RNA > nf.low & percent.mt < mt.cutoff & nFeature_RNA < nf.high)
}
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all.list <- SplitObject(all, split.by = group.by)

for (i in 1:length(all.list)){
  all.list[[i]] <- NormalizeData(all.list[[i]], verbose = FALSE)
  all.list[[i]] <- FindVariableFeatures(all.list[[i]], selection.method = "vst",nfeatures = nfeatures, verbose = FALSE)}

reference.list <- all.list
all.anchors <- FindIntegrationAnchors(object.list = reference.list, dims =1:dim.use)
all.integrated <- IntegrateData(anchorset = all.anchors, dims = 1:dim.use)
DefaultAssay(all.integrated) <- "integrated"
all.integrated <- ScaleData(all.integrated, verbose = FALSE)
npcs <- dim.use+10
all.integrated <- RunPCA(all.integrated, npcs = npcs, verbose = FALSE)

all.integrated <- FindNeighbors(all.integrated, reduction = "pca", dims = 1:dim.use)
all.integrated <- FindClusters(all.integrated, resolution = res)
all.integrated <- RunUMAP(all.integrated, reduction = "pca", dims = 1:dim.use)

return(all.integrated)
}

sct_integ <- function(obj,group.by=NULL,
                         mt.pattern="^MT-",mt.list=NULL,dim.use=30,mt.cutoff=5,
                         nf.low=500,nf.high=6000,nfeatures=3000,
                         res=1.5,is.mt=NULL,is.subset=NULL,is.pt=NULL, conserve.memory = TRUE) {
all <- obj
if (is.null(is.mt)) {
if (is.null(mt.list)) {
all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = mt.pattern)
}else{
mt.list <- mt.list[which(mt.list %in% rownames(all))]
all[["percent.mt"]] <- PercentageFeatureSet(all, features = mt.list)
}
}

if (!is.null(is.pt)) {
all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = "^pesudo-gene")
}


if (is.null(is.subset)) {
all <- subset(all, subset = nFeature_RNA > nf.low & percent.mt < mt.cutoff & nFeature_RNA < nf.high)
}

all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- all

ifnb.list <- SplitObject(obj, split.by = group.by)

if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform, vars.to.regress = "percent.mt", verbose = FALSE,conserve.memory = conserve.memory)
}else{
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform, vars.to.regress = "percent.mt", verbose = FALSE, method = "glmGamPoi",conserve.memory = conserve.memory)
}

features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = nfeatures)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
    anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:dim.use)

all.integrated <- immune.combined.sct
all.integrated <- FindNeighbors(all.integrated, reduction = "pca", dims = 1:dim.use)
all.integrated <- FindClusters(all.integrated, resolution = res)
all.integrated <- RunUMAP(all.integrated, reduction = "pca", dims = 1:dim.use)

return(all.integrated)
}

harmony_integ <- function(obj,group.by=NULL,
                         mt.pattern="^MT-",mt.list=NULL,dim.use=20,mt.cutoff=5,
                         nf.low=500,nf.high=6000,nfeatures=3000,
                         res=1.5,is.mt=NULL,is.subset=NULL) {
all <- obj
if (is.null(is.mt)) {
if (is.null(mt.list)) {
all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = mt.pattern)
}else{
mt.list <- mt.list[which(mt.list %in% rownames(all))]
all[["percent.mt"]] <- PercentageFeatureSet(all, features = mt.list)
}
}

if (is.null(is.subset)) {
all <- subset(all, subset = nFeature_RNA > nf.low & percent.mt < mt.cutoff & nFeature_RNA < nf.high)
}

all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(all)
all <- ScaleData(all , features = all.genes, vars.to.regress = "nCount_RNA")
saveRDS(all,"regress.rds")
all <- RunPCA(all, features = VariableFeatures(object = all))
all <- RunHarmony(all, group.by , plot_convergence = F,dims.use = 1:dim.use)
Combine <- all
Combine = RunTSNE(Combine, reduction = "harmony", dims = 1:dim.use)
Combine = RunUMAP(Combine, reduction = "harmony", dims = 1:dim.use)
Combine = FindNeighbors(Combine, reduction = "harmony",dims = 1:dim.use)
Combine = FindClusters(Combine, resolution = res)
return(Combine)
}

dimplot_2groups <- function(obj, 
                            group1="seurat_clusters", group2="orig.ident",
                            split.by=NULL,
                            wid=1800,hei=1000,
                            raster=NULL, pt.size=0.8,label.size=5,prefix='fig') {
if (is.null(raster)) {
png(paste0(prefix,"_",group1,"_",group2,".png"),width = 1800, height = 1000)
print(DimPlot(obj, reduction = "umap", label = T, label.size = label.size,repel =T,group.by=group1,pt.size=pt.size)+scale_color_manual(values=cor) -
                   (DimPlot(obj, reduction = "umap", label = T, group.by = group2, label.size = label.size,repel =T,pt.size=pt.size)+scale_color_manual(values=cor)))
dev.off()

if (!is.null(split.by)) {
png(paste0(prefix,"_",split.by, "_split.png"),width = wid, height = hei)
print(DimPlot(obj, reduction = "umap", label = T, split.by = split.by, label.size = label.size,repel =T, pt.size=pt.size)+scale_color_manual(values=cor))
dev.off()

}
}else{
png(paste0(prefix,"_",group1,"_",group2,".png"),width = 1800, height = 1000)
print(DimPlot(obj, reduction = "umap", label = T, label.size = label.size,repel =T,raster=FALSE,group.by=group1, pt.size=pt.size)+scale_color_manual(values=cor) -
                   (DimPlot(obj, reduction = "umap", label = T, group.by = group2, label.size = label.size,repel =T,raster=FALSE, pt.size=pt.size)+scale_color_manual(values=cor)))
dev.off()

if (!is.null(split.by)) {
png(paste0(prefix,"_",split.by, "_split.png"),width = wid, height = hei)
print(DimPlot(obj, reduction = "umap", label = T, split.by = split.by, label.size = label.size, repel =T,raster=FALSE, pt.size=pt.size)+scale_color_manual(values=cor))
dev.off()
}
}
}


re_cluster <- function(obj,res=1, group1="seurat_clusters",group2="orig.ident",split.by=NULL,wid=1800,hei=1000,assays=NULL) {
if (!is.null(assays)) {
DefaultAssay(obj) <- assays
}
obj <- FindClusters(obj, resolution = res)
dimplot_2groups(obj,group1=group1,group2=group2,split.by=NULL,wid=wid,hei=hei)
return(obj@meta.data)
}

top_cell  <- function(obj,top=10000) {
meta <- obj@meta.data
ncell <- dim(meta)[1]
if (ncell > top) {
meta <- meta[order(meta$nCount_RNA,decreasing=T),]
cellist <- rownames(meta)[1:top]
obj <- subset(obj,cells=cellist)
return(obj)

}else{
return(obj)
}

}

top_obj <- function(obj,top=10000,group.by="NULL") {
all.list <- SplitObject(obj, split.by = group.by)
seob.list <- lapply(all.list,top_cell,top=top)
obj <- Reduce(merge,seob.list)
return(obj)
}

hl_plot <- function(obj, group.by=NULL,raster=NULL) {
n <- sort(unique(obj@meta.data[,group.by]))
plist<-list()
Idents(obj) <- group.by
if (is.null(raster)) {
for (i in 1:length(n)) {
plist[[i]] <- DimPlot(obj, cells.highlight = rownames(obj@meta.data[which(obj@meta.data[,group.by] == n[i]),])) +ggtitle(n[i])+theme(legend.position= "none")
}
}else{
for (i in 1:length(n)) {
plist[[i]] <- DimPlot(obj, cells.highlight = rownames(obj@meta.data[which(obj@meta.data[,group.by] == n[i]),]),raster=F) +ggtitle(n[i])+theme(legend.position= "none")
}
}

hei=ceiling(length(n)/4)*750
p <- plist
jpeg(paste0("highlight_dimplot_",group.by,"_split.jpg"),width=3000,height=hei)
print(cowplot::plot_grid(plotlist=p,ncol=4))
dev.off()

}


tiff_dimplot_2groups <- function(obj, 
                            group1="seurat_clusters", group2="paper",
                            split.by=NULL,
                            wid=1800,hei=1000,
                            raster=NULL, label=T,dpi=NA,pt.size=0.2) {
if (!is.na(dpi)) {
n <- ceiling(dpi/300)*2000
hei1 <- n
wei1 <- 2*n
}else{
hei1=2000
wei1=4000
}

legend1 <- guides(color = guide_legend(override.aes = list(size=1), ncol=1) )

if (is.null(raster)) {
tiff(paste0(group1,"_",group2,".tiff"),width = wei1, height = hei1,res=dpi)
p1 <- DimPlot(obj, reduction = "umap", label = T, label.size = 2,repel =T,group.by=group1,pt.size=pt.size)+NoAxes()+scale_color_manual(values=cor)+legend1 
p2 <- DimPlot(obj, reduction = "umap", label = F, group.by = group1,pt.size=pt.size)+NoAxes()+scale_color_manual(values=cor)+legend1
p3 <- p1-p2
print(p3)

dev.off()

if (!is.null(split.by)) {
tiff(paste0(split.by, "_split.tiff"),width = wid, height = hei,res=dpi)
print(DimPlot(obj, reduction = "umap", label = T, split.by = split.by, label.size = 5,repel =T)+scale_color_manual(values=cor))
dev.off()
}
}else{
tiff(paste0(group1,"_",group2,".tiff"),width = wei1, height = hei1,res=dpi)
p1 <- DimPlot(obj, reduction = "umap", label = T, label.size = 2,repel =T,raster=FALSE,group.by=group1,pt.size=pt.size)+NoAxes+scale_color_manual(values=cor)+legend1
p2 <- DimPlot(obj, reduction = "umap", label = F, group.by = group1, raster=FALSE,pt.size=pt.size)+NoAxes()+scale_color_manual(values=cor)+legend1
p3 <- p1-p2
print(p3)

dev.off()

if (!is.null(split.by)) {
tiff(paste0(split.by, "_split.png"),width = wid, height = hei,res=dpi)
print(DimPlot(obj, reduction = "umap", label = T, split.by = split.by, label.size = 5,repel =T,raster=FALSE)+scale_color_manual(values=cor))
dev.off()
}
}
}

pdf_tiff_dimplot_2groups <- function(obj, 
                            group1="seurat_clusters", group2="paper",
                            split.by=NULL,
                            wid=1800,hei=1000,
                            raster=NULL, label=T,dpi=NA,pt.size=0.2) {
if (!is.na(dpi)) {
n <- ceiling(dpi/300)*2000
hei1 <- n
wei1 <- 2*n
}else{
hei1=2000
wei1=4000
}

legend1 <- NoLegend()

if (is.null(raster)) {
tiff(paste0(group1,"_",group2,".tiff"),width = wei1, height = hei1,res=dpi)
p1 <- DimPlot(obj, reduction = "umap", label = T, label.size = 2,repel =T,group.by=group1,pt.size=pt.size)+NoAxes()+scale_color_manual(values=cor)+legend1 
p2 <- DimPlot(obj, reduction = "umap", label = F, group.by = group1,pt.size=pt.size)+NoAxes()+scale_color_manual(values=cor)+legend1
p3 <- p1-p2
print(p3)

dev.off()
pdf(paste0(group1,"_",group2,".pdf"),width = wei1/100, height = hei1/100)
print(p3)
dev.off()


if (!is.null(split.by)) {
tiff(paste0(split.by, "_split.tiff"),width = wid, height = hei,res=dpi)
print(DimPlot(obj, reduction = "umap", label = T, split.by = split.by, label.size = 5,repel =T)+scale_color_manual(values=cor))
dev.off()
}
}else{
tiff(paste0(group1,"_",group2,".tiff"),width = wei1, height = hei1,res=dpi)
p1 <- DimPlot(obj, reduction = "umap", label = T, label.size = 2,repel =T,raster=FALSE,group.by=group1,pt.size=pt.size)+NoAxes+scale_color_manual(values=cor)+legend1
p2 <- DimPlot(obj, reduction = "umap", label = F, group.by = group1, raster=FALSE,pt.size=pt.size)+NoAxes()+scale_color_manual(values=cor)+legend1
p3 <- p1-p2
print(p3)

dev.off()

pdf(paste0(group1,"_",group2,".pdf"),width = wei1/100, height = hei1/100)
print(p3)
dev.off()

if (!is.null(split.by)) {
tiff(paste0(split.by, "_split.png"),width = wid, height = hei,res=dpi)
print(DimPlot(obj, reduction = "umap", label = T, split.by = split.by, label.size = 5,repel =T,raster=FALSE)+scale_color_manual(values=cor))
dev.off()
}
}
}

