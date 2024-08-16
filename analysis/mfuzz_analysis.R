library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 100000 * 1024^5)
library(Seurat)
library(Mfuzz)

plot_mfuzz <- function(obj,assay='RNA',mode='knn',nclusters=NULL,prefix='out',gene.use=NULL,nfeatures=3000,group.by='ident',inlevels=NULL){
Idents(obj) <- group.by
DefaultAssay(obj) <- assay
obj <- NormalizeData(obj)
if (is.null(inlevels)) {
inlevels <- sort(levels(obj))
}
obj@meta.data[,group.by] <- factor(obj@meta.data[,group.by],levels=inlevels)

if (is.null(gene.use)) {
obj <- FindVariableFeatures(obj,nfeatures=nfeatures)
gene.use <- VariableFeatures(obj)
}else{
c1 <- which(gene.use %in% rownames(obj))
gene.use <- gene.use[c1]
}
ave <- AverageExpression(obj,group.by = group.by,assays=c(assay),features=gene.use)
ingene <- ave$RNA
gene_tpm <- data.matrix(ingene)
eset <- new("ExpressionSet",exprs = gene_tpm)
gene.r <- filter.NA(eset, thres=0.25)
gene.f <- fill.NA(gene.r,mode=mode)
tmp <- filter.std(gene.f,min.std=0)
gene.s <- standardise(tmp)

if (is.null(nclusters)) {
nclusters=length(colnames(gene_tpm))
}
print(paste0('nclusters is:',nclusters))
c <- nclusters
m <- mestimate(gene.s)
cl <- mfuzz(gene.s, c = c, m = m)
saveRDS(cl,paste0(prefix,'.rds'))
write.table(cl$cluster,paste0(prefix,"_output.txt"),quote=F,row.names=T,col.names=F,sep="\t")

nrow <- ceiling(nclusters/3)

pdf(paste0(prefix,'.pdf'),25,ceiling(nclusters/6)*10)
mfuzz.plot(gene.s,cl,mfrow=c(nrow,3),new.window= FALSE)
mfuzz.plot2(gene.s,cl,mfrow=c(nrow,3),centre=TRUE,centre.lwd=4,time.labels=colnames(gene_tpm),x11 = FALSE)
dev.off()
}

c1<- which(all$first %in% c('hs','ma','mm'))
obj=subset(all,cells=Cells(all)[c1])

gene.use <- readLines('allTFs_hg38.txt')

prefix='layer'
group.by='third'
plot_mfuzz(obj,assay='RNA',mode='knn',nclusters=NULL,prefix=prefix,gene.use=gene.use,nfeatures=3000,group.by=group.by)

hvg=readRDS('11.ex_plot/hvg2000.rds')
plot_mfuzz(obj,assay='RNA',mode='knn',nclusters=NULL,prefix='hvg',gene.use=hvg,nfeatures=3000,group.by=group.by)
plot_mfuzz(obj,assay='RNA',mode='knn',nclusters=6,prefix='hvg_6cluster',gene.use=hvg,nfeatures=3000,group.by=group.by)
