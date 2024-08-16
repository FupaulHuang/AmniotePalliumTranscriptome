library(DoubletFinder)
library(Seurat)
library(dplyr)
library(RColorBrewer)
blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)

rm_doublet <- function(name=NULL,input=NULL,dim.usage=30,auto="true",gene.list=NULL) {
inpath <- list.files(path = input,pattern = name,full.names = T)
if (auto=="true") {
inpath <- paste0(inpath,"/","04.Matrix/")
}
EC <- Read10X(inpath,gene.column=1)
if (!is.null(gene.list)) {
genes <- rownames(EC)
genes <- intersect(gene.list,genes)
EC <- EC[genes,]
}
EC <- CreateSeuratObject(EC)
EC <- NormalizeData(EC)
EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = 3000)
EC <- ScaleData(EC)
EC <- RunPCA(EC)
EC <- RunUMAP(EC, dims = 1:dim.usage)

Find_doublet <- function(data){
sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
nExp_poi <- round(0.05*ncol(data))
p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
#data<-subset(data,subset=doublet_info=="Singlet")
return(data)
}

EC<-Find_doublet(EC)
EC<-subset(EC,subset=doublet_info=="Singlet")
EC@meta.data$library = name
c <- grep("pANN_",colnames(EC@meta.data))
EC@meta.data <- EC@meta.data[,-c]
print(paste0(name," cells: ", length(EC@meta.data$orig.ident)," is read!"," Time:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
return(EC)
}

seob_cluster <- function(obj,
                         mt.pattern="^MT-",mt.list=NULL,dim.use=30,mt.cutoff=10,
                         nf.low=200,nf.high=6000,nfeatures=3000,
                         res=1,is.mt=NULL,is.subset=T) {
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
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = nfeatures)
all.genes <- rownames(all)
all <- ScaleData(all, features = all.genes)

all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
all <- FindNeighbors(all, dims = 1:dim.use)
all <- FindClusters(all, resolution = res)
all <- RunUMAP(all, dims = 1:dim.use)
return(all)
}

plot_markers_split <- function(obj,gene.list=NULL,output="./",raster=NULL) {
all <- obj
gene <- unique(gene.list)
c <- which(gene %in% rownames(obj))
gene <- gene[c]
for (i in 1:length(gene)){
pf <- paste0(i,"_",gene[i],"_Featureplot.png")
pf <- gsub("/","_",pf)
png(paste0(output,pf))
if (is.null(raster)) {
print(FeaturePlot(obj,features = gene[i]))
}else{
print(FeaturePlot(obj,features = gene[i], raster=F))
}
dev.off()
}

}

subset_deg <- function(obj,group.by="seurat_clusters",sp.size=NULL,output="./",
                       min.pct=0.25,logfc.threshold=0.25,only.pos=F,assays ="RNA",order=F) {
all <- DietSeurat(obj)
if (!is.null(sp.size)) {
seob_list <- list()
i <- 1
for (sc in unique(all@meta.data[,group.by])){
cellist <- colnames(all)[which(all@meta.data[,group.by] == sc)]

ob <- subset(all, cells=cellist)

if (length(colnames(ob)) > sp.size) {
ob <- subset(ob,cells=sample(colnames(ob), sp.size))
}
seob_list[[i]] <- ob
i <- i+1
}

#n <- length(seob_list)
#all <- merge(seob_list[[1]], seob_list[c(2:n)])
all <- Reduce(merge,seob_list)
}

all_markers <- FindAllMarkers(all, only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold, verbose = F,assays = assays,order=order)
if (length(grep("avg_log2FC",colnames(all_markers))!=0)) {
all_markers$FC <- 2^(all_markers$avg_log2FC)
}else{
all_markers$FC <- exp(all_markers$avg_logFC)
}


write.table(all_markers, paste0(output, "deg_sample",sp.size,".xls"),sep="\t",quote=F)
}

library(cowplot)
library(RColorBrewer)
library(ggpubr)
library(ggsci)

colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))
len <- 100
cor<-c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"),colpalettes,rainbow(len))

sct_seob_cluster <- function(obj,
                         mt.pattern="^MT-",mt.list=NULL,dim.use=30,mt.cutoff=10,
                         nf.low=200,nf.high=6000,nfeatures=3000,
                         res=1,is.mt=NULL,is.subset=NULL,is.pt=NULL) {
pbmc <- obj
if (is.null(is.mt)) {
if (is.null(mt.list)) {
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = mt.pattern)
}else{
mt.list <- mt.list[which(mt.list %in% rownames(pbmc))]
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features = mt.list)
}
}

if (!is.null(is.pt)) {
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = '^pesudo_gene-')

}


if (is.null(is.subset)) {
pbmc <- subset(pbmc, subset = nFeature_RNA > nf.low & percent.mt < mt.cutoff & nFeature_RNA < nf.high)
}
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = nfeatures)
if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
}else{
pbmc <- SCTransform(pbmc, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
}

pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 40, verbose = F)
pbmc <- RunUMAP(pbmc, dims = 1:dim.use, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:dim.use, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = res)
return(pbmc)
}

rm_doublet_top <- function(name=NULL,input=NULL,dim.usage=30,auto="true",umi=NULL,topumi=10000,gene.list=NULL) {
inpath <- list.files(path = input,pattern = name,full.names = T)
if (auto=="true") {
inpath <- paste0(inpath,"/","04.Matrix/")
}
EC <- Read10X(inpath,gene.column=1)
if (!is.null(gene.list)) {
genes <- rownames(EC)
genes <- intersect(gene.list,genes)
EC <- EC[genes,]
}
EC <- CreateSeuratObject(EC)
if (!is.null(umi)){
EC <- subset(EC,subset = nCount_RNA > umi)
}
if (dim(EC)[1] > topumi) {
meta <- EC@meta.data
meta <- meta[order(meta$nCount_RNA,decreasing=T),]
meta <- meta[1:topumi,]
cellist <- rownames(meta)
EC <- subset(EC,cells=cellist)
}
EC <- NormalizeData(EC)
EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = 3000)
EC <- ScaleData(EC)
EC <- RunPCA(EC)
EC <- RunUMAP(EC, dims = 1:dim.usage)

Find_doublet <- function(data){
sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
nExp_poi <- round(0.05*ncol(data))
p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
return(data)
}
EC<-Find_doublet(EC)
EC<-subset(EC,subset=doublet_info=="Singlet")
EC@meta.data$library = name
c <- grep("pANN_",colnames(EC@meta.data))
EC@meta.data <- EC@meta.data[,-c]
print(paste0(name," cells: ", length(EC@meta.data$orig.ident)," is read!"," Time:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
return(EC)
}

qc_vlnplot <- function(obj, assay =NULL,
                       features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                       group.by="seurat_clusters",pt.size=0,
                       pf=NULL,wid=1800,hei=1000,raster=NULL,stack=TRUE,flip=TRUE){
p1 <- VlnPlot(obj, features = features, group.by = group.by, pt.size = pt.size,assay=assay, stack=stack,flip=flip)
if (is.null(pf)) {
pf <- paste0("vlnplot_",group.by,".jpg")
}
jpeg(pf,wid,hei)
print(p1)
dev.off()

if (is.null(raster)) {
p1 <- FeaturePlot(obj,features = features,ncol=3)
}else{
p1 <- FeaturePlot(obj,features = features, raster=F,ncol=3)
}
if (is.null(pf)) {
pf <- paste0("featureplot_",group.by,".jpg")
}else{
pf <- paste0("featureplot_",pf)
}
jpeg(pf,2100,1000)
print(p1)
dev.off()
}


sample_seob <- function(obj,group.by="seurat_clusters",sp.size=NULL,diet="true",sp.total=1000) {
all <- obj
if (diet=="true") {
all <- DietSeurat(all,dimreducs = c('pca','umap'))
}

if (is.null(sp.size)) {
nlen <- length(unique(all@meta.data[,group.by]))
sp.size <- ceiling(sp.total/nlen)

}
seob_list <- list()
i <- 1
for (sc in unique(all@meta.data[,group.by])){
cellist <- colnames(all)[which(all@meta.data[,group.by] == sc)]
ob <- subset(all, cells=cellist)
if (length(colnames(ob)) > sp.size) {
ob <- subset(ob,cells=sample(colnames(ob), sp.size))
}
seob_list[[i]] <- ob
i <- i+1
}
all <- Reduce(merge,seob_list)

return(all)
}

Sample_seob <- function(obj,group.by="seurat_clusters",sp.size=NULL,diet="true",sp.total=1000) {
all <- obj
if (diet=="true") {
all <- DietSeurat(all,dimreducs = c('pca','umap'))
}

if (is.null(sp.size)) {
nlen <- length(unique(all@meta.data[,group.by]))
sp.size <- ceiling(sp.total/nlen)

}
ncellist <- c()
for (sc in unique(all@meta.data[,group.by])){
cellist <- colnames(all)[which(all@meta.data[,group.by] == sc)]
if (length(cellist) > sp.size) {
cellist=sample(cellist, sp.size)
}
ncellist <- c(ncellist,cellist)
}
all <- subset(all,cells=ncellist)
return(all)
}


markers_plot <- function(obj, assay ="RNA",
                       features = features,
                       group.by="seurat_clusters",pt.size=0,
                       pf=NULL,wid=1800,hei=1000,raster=NULL){
p1 <- VlnPlot(obj, features = features, group.by = group.by, pt.size = 0,assay=assay, stack=TRUE,flip=TRUE)
if (is.null(pf)) {
pf <- paste0("markers_vlnplot_",group.by,".jpg")
}
ngene <- length(features)
hei <- ceiling(ngene/10)*1000
jpeg(pf,wid,hei)
print(p1)
dev.off()


nrow=ceiling(ngene/3)
if (is.null(raster)) {
p1 <- FeaturePlot(obj,features = features,ncol=3, pt.size=pt.size)
}else{
p1 <- FeaturePlot(obj,features = features, raster=F,ncol=3,pt.size=pt.size)
}
if (is.null(pf)) {
pf <- paste0("markers_featureplot_",group.by,".jpg")
}else{
pf <- paste0("markers_featureplot_",pf)
}
hei <- 700*nrow
jpeg(pf,2100,hei)
print(p1)
dev.off()

p1 <- DotPlot(obj, features = features, group.by = group.by, assay=assay)
pf <- paste0("markers_dotplot_",group.by,".jpg")
ncluster <- length(unique(obj$seurat_clusters))
hei <- ceiling(ncluster/40)*1000
jpeg(pf,wid,hei)
print(p1)
dev.off()

}

doplot_deg <- function(obj,output="./",group.by="seurat_clusters",min.pct=0.25,logfc.threshold=0.25,only.pos=T,
                       assays='RNA',order=F,ntop=3) {
#tmp <- obj
obj <- DietSeurat(obj)
Idents(obj) <- group.by
obj <- subset(obj,downsample = 5000)

all_markers <- FindAllMarkers(obj, only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold, verbose = F,assays = assays,order=order)
#nwt <- colnames(all_markers)[grep("avg",colnames(all_markers))]

all_markers %>% group_by(cluster) %>% top_n(ntop,wt=avg_log2FC) %>% data.frame() -> df1
gene <- unique(df1$gene) 

p1 <- DotPlot(obj, features = gene, group.by = group.by, assay=assays)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
pf <- paste0("deg_dotplot_",group.by,".jpg")
ncluster <- length(unique(obj@meta.data[,group.by]))
hei <- ceiling(ncluster/40)*1200
wid <- ceiling(length(gene)/40)*1800
jpeg(pf,wid,hei)
print(p1)
dev.off()
}

pct_data <- function(meta,x,fill){
meta$g1 <- meta[,x]
meta$g2 <- meta[,fill]
meta %>%  group_by(g1) %>%  count(g2) %>% data.frame() -> df2
v1 <- table(meta$g1)
v2 <- table(df2$g1)

total <- c()
for (n in names(v1)) {
    total <- c(total,rep(v1[n],v2[n]))
}

df2$total <- total
df2$pct <- df2$n/df2$total
write.table(df2,paste0(x,"_",fill,"_pct.xls"),sep='\t',quote=F)
}

pct_barplot <- function(meta,x,fill,wid=10,hei=10){
meta$tech3 <- meta[,x]
meta$Celltype <- meta[,fill]

pdf(paste0(x,"_",fill,"_fraction_compare.pdf"),width = wid, height = hei)
pA<-ggplot(meta,aes(x=tech3,fill=Celltype))+geom_bar(stat="count",position="fill")+ theme_bw()+coord_flip()+scale_fill_manual(values = cor)+
    theme(axis.title.x = element_text(size = 20,face = "bold"),
    axis.title.y = element_text(size = 10,face = "bold")) +
    theme(axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 8))
print(pA)
dev.off()
pct_data(meta=meta,x=x,fill=fill)
}

diet_obj <- function(obj,assays=NULL,dims=NULL,sct='true',is.save=T){
if (is.null(assays)){
assays <- Assays(obj)
}
if (is.null(dims)){
dims <- Reductions(obj)
}
if (sct=='true'){
sct <- obj@assays$SCT
assays <- assays[-grep("SCT",assays)]
if (is.save) {
saveRDS(sct,'sct.rds')
}
}
obj <- DietSeurat(obj, assays=assays,dimreducs=dims)
if (is.save) {
saveRDS(obj,'diet.rds')
meta <- obj@meta.data
write.table(meta, paste0("metadata.xls"), sep="\t", quote=F)
}
return(obj)
}


