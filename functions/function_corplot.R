
library(ggplot2)
library(Seurat)
library(corrplot)
library(stringr)

cortable <- function(ref=list(c(pbmc_small,"RNA","groups","label1")),
                    que=list(c(pbmc_small,"RNA","groups","label2")),
                    features=NULL,overlap="false",is.rowmeans=NULL,corr.method='spearman') {

###FUN1 AVERAGEEXPRESSION
preave <- function(obj,
                   features=c(rownames(pbmc_small)),
                   group="groups",
                   assay="RNA",
                   label="label1") {
DefaultAssay(obj) <- assay
ob1 <- NormalizeData(obj)
grp1 <- group
Idents(ob1)<-ob1[[grp1]]
ave1 <- AverageExpression(ob1,features = features,assays=assay)
ave1 <- ave1[[1]]
colnames(ave1) <- paste0(label,"_",colnames(ave1))

if (is.null(is.rowmeans)) {
Sp1 = ave1
  avg = rowMeans(Sp1)
  Sp1 = sweep(Sp1,1,avg,"/")
  rm(avg)
}else{
Sp1 = ave1
}

Sp1[is.nan(Sp1)] <- 0
ave1 <- Sp1
return(ave1)
}
###FUN2 GET FEATURES
getfeatures <- function(ref=list(c(pbmc_small,"RNA","groups","label1")),
                        que=list(c(pbmc_small,"RNA","groups","label2"))) {
lenref <- length(ref)
lenque <- length(que)

reflist <- list()
for (i in 1:lenref) {
tmp <- ref[[i]][[1]]
DefaultAssay(tmp) <- ref[[i]][[2]]
reflist[[i]] <- tmp
}
geneset1 <- lapply(reflist[1:lenref],rownames)
gene1 <- Reduce(intersect, geneset1)

quelist <- list()
for (i in 1:lenque) {
tmp <- que[[i]][[1]]
DefaultAssay(tmp) <- que[[i]][[2]]

quelist[[i]] <- tmp
}
geneset2 <- lapply(quelist[1:lenque],rownames)
gene2 <- Reduce(intersect, geneset2)

final_featues <- intersect(gene1,gene2)
return(final_featues)
}
if (is.null(features)) {
features <- getfeatures(ref,que)
}else{
features1 <- getfeatures(ref,que)
features <- intersect(features,features1)
}
write.table(features,"use_gene.list",sep='\t',quote=F,row.names=F,col.names=F)


###FUN3 MERGE AVERAGEEXPRESSION
preave_list <- function(inlist, features=NULL) {
len <- length(inlist)
matlist <- list()
for (i in 1:length(inlist)) {
matlist1 <- preave(obj=inlist[[i]][[1]],
                  features=features,
                  assay=inlist[[i]][[2]],                       
                  group=inlist[[i]][[3]],                        
                  label=inlist[[i]][[4]])                        
matlist[[i]] <- matlist1
}
mat <- Reduce(cbind,matlist)
return(mat)
#lapply(inlist,preave(),features=features,group=)
}

###Compute cor values
Sp1 = preave_list(ref,features=features)
colnames(Sp1) <- paste0(colnames(Sp1),"_ref")
Sp2 = preave_list(que,features=features)
colnames(Sp2) <- paste0(colnames(Sp2),"_que")

geTable = merge(Sp1,Sp2, by='row.names', all=F)

rownames(geTable) = geTable$Row.names
geTable = geTable[,2:ncol(geTable)]
# corr.method = c('spearman', 'pearson') etc.
Corr.Coeff.Table = cor(geTable,method=corr.method)

###Estimate p-value
nPermutations = 1000
shuffled.cor.list = list()
  pb   <- txtProgressBar(1, 100, style=3)

  for (i in 1:nPermutations){
    shuffled = apply(geTable[,1:ncol(Sp1)],1,sample)
    shuffled2 = apply(geTable[,(ncol(Sp1)+1):ncol(geTable)],1,sample)
    shuffled = cbind(t(shuffled),t(shuffled2))
    shuffled.cor = cor(shuffled,method=corr.method)
    shuffled.cor.list[[i]] = shuffled.cor
    rm(list=c('shuffled','shuffled2','shuffled.cor'))
    if ((i %% 100) ==0){
      setTxtProgressBar(pb, (i*100)/nPermutations)
    }
  }

  p.value.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))

  rownames(p.value.table) = colnames(geTable)
  colnames(p.value.table) = colnames(geTable)

  shuffled.mean.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
  rownames(shuffled.mean.table) = colnames(geTable)
  colnames(shuffled.mean.table) = colnames(geTable)

  a = combn(1:ncol(geTable),2)
  for (i in 1:ncol(a)){
    cor.scores = sapply(shuffled.cor.list,"[",a[1,i],a[2,i])
    shuffled.mean.table[a[1,i],a[2,i]] = mean(cor.scores)
    shuffled.mean.table[a[2,i],a[1,i]] = mean(cor.scores)
    p.value = mean(abs(cor.scores)>=abs(Corr.Coeff.Table[a[1,i],a[2,i]]))
    p.value.table[a[1,i],a[2,i]] = p.value
    p.value.table[a[2,i],a[1,i]] = p.value
    rm(list=c('cor.scores','p.value'))
    setTxtProgressBar(pb, (i*100)/ncol(a))
  }
p.value.table[is.na(p.value.table)] <- 1
if (overlap=="false") {
M <- p.value.table
mat <- M[,grep('_ref',colnames(M))]
mat <- mat[grep('_que',rownames(M)),]
p.value.table  <- mat

M <- Corr.Coeff.Table
mat <- M[,grep('_ref',colnames(M))]
mat <- mat[grep('_que',rownames(M)),]
Corr.Coeff.Table <- mat
}

return(list(Corr.Coeff.Table,p.value.table))
}


corplot <- function(cor.table=NULL,
                    pva.table=NULL,
                    cutoff=0,
                    pf=NULL,
                    col=colorRampPalette(c("darkblue", "white","darkred")),
                    order="original",
                    wid=1000,
                    hei=1000,
                    label.size=1) {

cor.table[cor.table<cutoff] <- 0

if (is.null(pf)) {
pf <- paste0("corrplot_filter",cutoff,"_",order,".jpg")
}
jpeg(pf,wid,hei)
print(
      corrplot(cor.table, order=order,tl.pos="lt", method="color", tl.col="black",cl.lim=c(min(cor.table),max(cor.table)), is.corr=F,tl.cex=label.size, sig.level=(-0.05),insig="pch", pch=19, pch.cex=0.25,pch.col="black",p.mat=pva.table, col=col(200), main= paste(pf),mar=c(3,1,5,1),cl.align.text="l")
      )
dev.off()

}

library(pheatmap)
cor_heatmap <- function(cor.table,pva.table=NULL,
                        col= colorRampPalette(c("darkblue", "white","darkred"))(256),
                        pf=NULL,wid=1000,hei=1000,res=1200,
                        cutoff=-1,
                        scale="none",is.pdf=NULL){

if (cutoff!=-1){
cor.table[cor.table<cutoff] <- 0
}
if (is.null(pva.table)) {
p1<-pheatmap(cor.table,scale=scale,show_colnames = T,show_rownames = T,fontsize = 40,
             cluster_rows = T,cluster_cols = T,border_color = "NA",color = col,res=res,filename=NA)
}else{
pmt <- pva.table
ssmt <- pmt< 0.01
pmt[ssmt] <-'**'
smt <- pmt >0.01& pmt <0.05
pmt[smt] <- '*'
pmt[!ssmt&!smt]<- ''

p1<-pheatmap(cor.table,scale=scale,show_colnames = T,show_rownames = T,fontsize = 40,
             cluster_rows = T,cluster_cols = T,border_color = "NA",color = col, res=res,filename=NA,
             display_numbers = pmt, number_color = "black")

}
if (is.null(pf)) {
pf <- paste0("corheatmap_filter",cutoff,".jpg")
}
jpeg(pf,wid,hei)
print(p1)
dev.off()

if (!is.null(is.pdf)) {
pf <- paste0("corheatmap_filter",cutoff,".pdf")
pdf(pf,wid/100,hei/100)
print(p1)
dev.off()
}

}

  sh_cor <- function(i,Sp1,geTable) {
    shuffled = apply(geTable[,1:ncol(Sp1)],1,sample)
    shuffled2 = apply(geTable[,(ncol(Sp1)+1):ncol(geTable)],1,sample)
    shuffled = cbind(t(shuffled),t(shuffled2))
    shuffled.cor = cor(shuffled,method=corr.method)
    shuffled.cor.list[[i]] = shuffled.cor
    return(shuffled.cor)
    rm(list=c('shuffled','shuffled2','shuffled.cor'))
    if ((i %% 100) ==0){
      setTxtProgressBar(pb, (i*100)/nPermutations)
    }
  }
  #shuffled.cor.list <- lapply(c(1:nPermutations),sh_cor,Sp1=Sp1,geTable=geTable)

# https://github.com/satijalab/seurat/issues/2617
# RenameGenesSeurat  ------------------------------------------------------------------------------------

RenameGenesSeurat <- function(obj,newnames,gene.use=NULL,de.assay="Spatial") { 
  # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. 
  # It only changes obj@assays$RNA@counts, @data and @scale.data.
print("Run this before integration. It only changes obj@assays$*@counts, @data and @scale.data, @var.features,@reductions$pca@feature.loadings")
lassays <- Assays(obj)
#names(obj@assays)
try(assay.use <- obj@reductions$pca@assay.used)
if (!exists('assay.use')) { assay.use=de.assay}
DefaultAssay(obj) <- de.assay
if (is.null(gene.use)) {
all_genenames <- rownames(obj)
}else{
all_genenames <- gene.use
obj <- subset(obj,features=gene.use)
}

order_name <- function(v1,v2,ref){
v2 <- make.names(v2,unique=T)
df1 <- data.frame(v1,v2)
rownames(df1) <- df1$v1
df1 <- df1[ref,]
return(df1)
}

df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
all_genenames <- df1$v1
newnames <- df1$v2

if ('SCT' %in% lassays) {
if ('SCTModel.list' %in%  slotNames(obj@assays$SCT)) {
obj@assays$SCT@SCTModel.list$model1@feature.attributes <- obj@assays$SCT@SCTModel.list$model1@feature.attributes[all_genenames,]
rownames(obj@assays$SCT@SCTModel.list$model1@feature.attributes) <- newnames
}
}

change_assay <- function(a1=de.assay,obj,newnames=NULL,all_genenames=NULL){
RNA <- obj@assays[a1][[1]]  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@var.features)) {
        df1 <- order_name(v1=all_genenames,v2=newnames,ref=RNA@var.features)
        all_genenames1 <- df1$v1
        newnames1 <- df1$v2
        RNA@var.features            <- newnames1
    }
    if (length(RNA@scale.data)){
        df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(RNA@scale.data))
        all_genenames1 <- df1$v1
        newnames1 <- df1$v2
        rownames(RNA@scale.data)    <- newnames1
    }

  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays[a1][[1]] <- RNA
  return(obj)
}

for (a in lassays) {
DefaultAssay(obj) <- a
df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
all_genenames1 <- df1$v1
newnames1 <- df1$v2
obj <- change_assay(obj=obj,a1=a,newnames=newnames1,all_genenames=all_genenames1)
}

hvg <- VariableFeatures(obj,assay=assay.use) 
if (length(obj@reductions$pca)){
    df1 <- order_name(v1=all_genenames,v2=newnames,ref=hvg)
    df1 <- df1[rownames(obj@reductions$pca@feature.loadings),]
    all_genenames1 <- df1$v1
    newnames1 <- df1$v2
    rownames(obj@reductions$pca@feature.loadings) <- newnames1
}
try(obj[[de.assay]]@meta.features <- data.frame(row.names = rownames(obj[[de.assay]])))
return(obj)
}



run_corplot <- function(ob1,ob3,ig1=NULL,ig2=NULL,col1=NULL,col2=NULL,lab1=NULL,lab2=NULL,gene.list=NULL,assay1='Spatial',assay2='Spatial',deg=NULL,hvg=NULL,invert=NULL,assay=NULL){

DefaultAssay(ob1) <-assay1
DefaultAssay(ob3) <-assay2   
if (!is.null(ig1)) {
ingene <- read.table(ig1,sep='\t',header=T,stringsAsFactor=F)
ob1 <- RenameGenesSeurat(obj=ob1,newnames= ingene[,1],gene.use=ingene[,3])
}

if (!is.null(ig2)) {
ingene <- read.table(ig2,sep='\t',header=T,stringsAsFactor=F)
ob3 <- RenameGenesSeurat(obj=ob3,newnames= ingene[,1],gene.use=ingene[,3])
}

if (!is.null(invert)) {
DefaultAssay(ob1)=assay1
DefaultAssay(ob3)=assay2
gene.list1=intersect(rownames(ob1), rownames(ob3))
c1 <- which(gene.list1 %in% gene.list)
len_list=length(c1)
print(paste0('len: ',len_list))
gene.list <- gene.list1[-c1]
}

que <- list(c(ob1,assay1,col1,lab1))
ref <- list(c(ob3,assay2,col2,lab2))

if (is.null(gene.list)) {
gene.list=intersect(rownames(ob1), rownames(ob3))
}

if (!is.null(deg)) {
Idents(ob1) <- col1
Idents(ob3) <- col2
deg1 <- FindAllMarkers(ob1,only.pos=F,min.pct = 0.25,logfc.threshold = 0.25,verbose = F)
deg2 <- FindAllMarkers(ob3,only.pos=F,min.pct = 0.25,logfc.threshold = 0.25,verbose = F)
gene.list1 <- unique(deg1$gene,deg2$gene)
gene.list <- intersect(gene.list1,gene.list)
}

if (!is.null(hvg)) {
Idents(ob1) <- col1
Idents(ob3) <- col2
ob1 <- FindVariableFeatures(ob1,nfeatures=hvg)
ob3 <- FindVariableFeatures(ob3,nfeatures=hvg)

gene.list1 <- unique(c(VariableFeatures(ob1),
                      VariableFeatures(ob3)))
#gene.list1 <- intersect(VariableFeatures(ob1),VariableFeatures(ob3))
c1=which(gene.list1 %in% gene.list)
gene.list1=gene.list1[c1]
if (!is.null(hvg)) {
gene.list1=gene.list1[1:len_list]
}
gene.list <- intersect(gene.list1,gene.list)
}


tmp <- cortable(ref,que,features=gene.list,corr.method='spearman')
saveRDS(tmp,"final_result.rds")
Corr.Coeff.Table <- tmp[[1]]
p.value.table <- tmp[[2]]

mat <- Corr.Coeff.Table
c1 <- mat>=0.999999
if (length(which(c1)) !=0) {
mat[c1] <- 0
imax <- max(mat)
mat[c1] <- imax
Corr.Coeff.Table <- mat
}

comp_table.species1 <- Corr.Coeff.Table
p_table.species1 <- p.value.table

corplot(cor.table=Corr.Coeff.Table,pva.table=p.value.table,cutoff=0,hei=1800)
corplot(cor.table=Corr.Coeff.Table,pva.table=p.value.table,cutoff=0.5,hei=1800)
cor_heatmap(Corr.Coeff.Table,p.value.table,wid=1500,hei=1500,is.pdf=T)
}

