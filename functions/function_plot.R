source("functions/function_merge.R")
source('functions//function_integ.R')
source('functions/function_heatmap.R')
library(stringr)
library(data.table)
library(tidyverse)
library(ggrepel)
library(RColorBrewer)

colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))
len <- 100
col<-c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"),colpalettes,rainbow(len))
cor <- col
blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)


dimplot_out <- function(obj,label='umap',corl=cor,pt.size=1.5,group1='sub',wei=1000,hei=1000){
obj$sub <- obj@meta.data[,group1]
jpeg(paste0(label,'_',group1,'.jpg'),wei,hei)
p1 <- DimPlot(obj, reduction = "umap", label = T, label.size = 0,repel =T,group.by=group1,pt.size=pt.size,raster=F)+scale_color_manual(values=corl)+NoAxes()+NoLegend() 
print(p1)
dev.off()

pdf(paste0(label,'_',group1,'_label','.pdf'), wei/100, hei/100)
p2 <- DimPlot(obj, reduction = "umap", label = F,group.by=group1, pt.size=0,raster=F)+scale_color_manual(values=corl)+NoAxes()+ggtitle(label)
print(p2)
dev.off()
}


plot_umap <- function(inrds,inmeta,label){
ob1 <- readRDS(inrds)
df1 <- read.table(inmeta,sep='\t',header=T)
ob1 <- subset(ob1,cells=rownames(df1))
ob1@meta.data <- df1
c1 <- which(ltype %in% unique(ob1$Sub))
ltype <- ltype[c1]
lcor <- lcor[ltype]
ob1$Sub <- factor(ob1$Sub,levels=ltype)

dimplot_out(ob1,label=label,corl=lcor,pt.size=1.5,group1='Sub',wei=1000,hei=1000)
}



pct_barplot <- function(meta,group.by,split.by,label) {
fn <- paste0(label,"_",group.by,"_pct.pdf")
pf <- fn 

meta1 <- meta
meta1$Celltype <- meta1[,group.by]
meta1$tech3 <- meta1[,split.by]
p1<-ggplot(meta1,aes(x=tech3,fill=Celltype))+geom_bar(stat="count",position="fill")+ theme_bw()+scale_fill_manual(values = cor)+
    xlab(label)+
    ylab("Ratio")+
    theme(axis.title.x = element_text(size = 20,face = "bold"),
    axis.title.y = element_text(size = 10,face = "bold")) +
    theme(axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 8))

pdf(pf,width = 18, height = 12)
print(p1)
dev.off()

for (m in unique(meta1$tech3)) {
    te <- subset(meta1,tech3==m)
    full <- te$Celltype
    len <- length(full)
    n <- 1
    ctn <- c()
    pctn <- c()
    countn <- c()
    for (i in unique(te$Celltype)) {
    count <- length(which(full %in% i))
    pct <- count/len
    ctn[n] <- i
    pctn[n] <- pct
    countn[n] <- count
    n <- n + 1
    }
    df <- data.frame(ctn,pctn,countn)
    df$tech <- m
    colnames(df) <- c("Celltype","pct","Count","individual")
    nam <- paste0("df_",m)
    assign(nam,df)
}
df <- get(paste0("df_",m))
for (i in unique(meta1$tech3)) {
    nam <- paste0("df_",i)
    df <- rbind(df,get(nam))
    df <- unique(df)
}

write.table(df,paste0(label,"_",group.by,"_fraction.xls"),sep="\t",quote=F,row.names=F)
}


qc_violin_facet <- function(meta,label='',group,lqc=c("nFeature_RNA","nCount_RNA" )){
fn <- paste0(label,"_violin_qc.pdf")
pf <- fn
meta$group <- meta[,group]

pdf(pf,5,5*length(unique(meta$group)))
for (i in lqc) {
meta$fake <- i
pb1 <- ggviolin(meta,x='fake',y=i,fill=group,add = "boxplot",width=1,add.params = list(width=0.1,size=1),facet.by = c(group))+
    xlab(i)+scale_y_log10()+
    scale_fill_manual(values=cor)+
    theme(axis.title.x = element_text(size = 20,face = "bold"),
    axis.title.y = element_text(size = 20,face = "bold")) +
    theme(axis.text.x = element_text(size = 15,),
    axis.text.y = element_text(size = 15,))

pb1 <- pb1 + facet_wrap(~group,ncol = 1) 

print(pb1)
}

dev.off()

nc <- c()
nf <- c()
nam <- c()
n <- 1
for (i in unique(meta[,group])){
    ob <- subset(meta,group==i)
    nc[n] <- median(ob$nCount_RNA)
    nf[n] <- median(ob$nFeature_RNA)
    nam[n] <- i
    n <- n + 1
}

df <- data.frame(nam,nc,nf)
colnames(df) <- c("Droup","UMI","Genes")
write.table(df,paste0(label,"_vlnplot_QC_median.xls"),sep="\t",quote=F)
}

qc_violin <- function(meta,label='',group,lqc=c("nFeature_RNA","nCount_RNA" )){
fn <- paste0(label,"_violin_qc.pdf")
pf <- fn
meta$group <- meta[,group]

meta$nFeature_RNA <- log10(meta$nFeature_RNA)
meta$nCount_RNA <- log10(meta$nCount_RNA)
pdf(pf,6,2*length(unique(meta$group)))
for (i in lqc) {
pb1 <- ggviolin(meta,x=group,y=i,fill=group,add = "boxplot",width=1,add.params = list(width=0.1,size=0.5),orientation = "horiz")+
    xlab(i)+
    #scale_y_log10()+
    scale_fill_manual(values=cor)+
    theme(axis.title.x = element_text(size = 20,face = "bold"),
    axis.title.y = element_text(size = 20,face = "bold")) +
    theme(axis.text.x = element_text(size = 15,),
    axis.text.y = element_text(size = 15,))

print(pb1)
}

dev.off()

nc <- c()
nf <- c()
nam <- c()
n <- 1
for (i in unique(meta[,group])){
    ob <- subset(meta,group==i)
    nc[n] <- median(ob$nCount_RNA)
    nf[n] <- median(ob$nFeature_RNA)
    nam[n] <- i
    n <- n + 1
}

df <- data.frame(nam,nc,nf)
colnames(df) <- c("Droup","UMI","Genes")
write.table(df,paste0(label,"_vlnplot_QC_median.xls"),sep="\t",quote=F)
}

pmarkers=function (obj, gene, col = c("lightgrey", "red"),pt.size=1) 
{
    for (i in 1:length(gene)) {
        png(paste0(i, "_", gene[i], "_Featureplot.png"))
        print(FeaturePlot(obj, features = gene[i], cols = col, pt.size=pt.size,
            raster = FALSE))
        dev.off()
        pdf(paste0(i, "_", gene[i], "_Featureplot.pdf"))
        print(FeaturePlot(obj, features = gene[i], cols = col, pt.size=pt.size,
            raster = FALSE))
        dev.off()
    }
}

plot_Dotplot=function (ob1, gene, ingroup = "seurat_clusters", inlabel = "out", 
    inassay = "RNA", mwid = 450, mhei = 450) 
{
    DefaultAssay(ob1) <- "RNA"
    gene = unique(gene)
    blu <- colorRampPalette(brewer.pal(9, "Blues"))(100)
    p1 = DotPlot(ob1, features = gene, group.by = ingroup) + 
        coord_flip() + scale_color_gradientn(colours = blu) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
            vjust = 0.5)) + ggtitle(inlabel)
    inwid = max(192 + length(unique(ob1@meta.data[, ingroup])) * 
        32, mwid)
    inlen = max(c(20 * length(gene), mhei))
    jpeg(paste0(inlabel, "_markers_dotplot.jpg"), inwid, inlen)
    print(p1)
    dev.off()
    pdf(paste0(inlabel, "_markers_dotplot.pdf"), inwid/100, inlen/100)
    print(p1)
    dev.off()
}

plot_vlnplot=function (ob1, gene, ingroup = "seurat_clusters", inlabel = "out", 
    inassay = "RNA", corl = c(brewer.pal(8, "Dark2"), brewer.pal(12, 
        "Paired"), brewer.pal(8, "Set2"), brewer.pal(9, "Set1"), 
        colpalettes, rainbow(20)),mwid = 500, mhei = 500) 
{
    DefaultAssay(ob1) <- "RNA"
    gene = unique(gene)
    p1 = VlnPlot(ob1, features = gene, group.by = ingroup, stack = T, 
        flip = T, fill.by = "ident") + scale_fill_manual(values = corl) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
            vjust = 0.5)) + NoLegend() + ggtitle(inlabel)
    inwid = max(length(unique(ob1@meta.data[, ingroup])) * 28, mwid)
    inlen = max(19 * length(gene), mhei)
    pdf(paste0(inlabel, "_markers_vlnplot.pdf"), inwid/100, inlen/100)
    print(p1)
    dev.off()
    jpeg(paste0(inlabel, "_markers_dotplot.jpg"), inwid, inlen)
    print(p1)
    dev.off()
}

deg_volcano_plot <- function(all_markers,label='out',fc.filter=1.5,p.val=0.01, wt='avg_log2FC',gene.col='gene',lb.filter=2){
try({
  a <- all_markers
  a$gene <- a[,gene.col]
  a$logfc <- a[,wt]
  if (wt=='avg_log2FC'){
  a$FC<-2^a$logfc
  }else{
  a$FC<-exp(a$logfc)
  }
  a$log2FC <- log2(a$FC)
  a_h<-a[a$FC>=fc.filter & a$p_val_adj<=p.val,]
  a_l<-a[a$FC<1/fc.filter & a$p_val_adj<p.val,]

  gene_list <- c(a_h$gene, a_l$gene)
  a_all<-rbind(a_h,a_l)
  pos<- a$gene %in% a_all$gene
  a_no<-a[!pos,]
  try(a_h$type<-"Up")
  try(a_l$type<-"Down")
  try(a_no$type<-"None")
  a_all<-rbind(a_h,a_no,a_l)
  mat<-a_all
  ran<-runif(nrow(a_all),min = 25,max =50)
  ran<-100^-ran
  mat$p_val_adj<-mat$p_val_adj+ran
  cdn <- which(mat$gene %in% gene_list)
  mat1 <- mat[cdn,]
pdf(paste0(label,"_vocalno_plot_",".pdf"),12,12)
  p1<-ggplot(mat,aes(log2(FC),-1*log10(p_val_adj)))+
            geom_point(aes(color =type))+
            geom_text_repel(data=subset(mat1, abs(mat1$log2FC) >=log2(lb.filter)),aes(x=log2(FC),y=-log10(p_val_adj),label=gene),size=5,point.padding = 0.5)+
            scale_color_manual(values=c('#377EB8',"lightgrey",'#E41A1C'))+
            theme(panel.background =element_blank(),axis.line=element_line(colour="black"))

  p1<-p1+theme(legend.title=element_blank(),legend.background = element_blank(),legend.key = element_blank(), legend.text = element_text(size=20), legend.position="NA")
  p1<-p1+theme_bw()+theme(axis.text = element_text(size = 20))
  p1<-p1+geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.8)+geom_hline(yintercept = -log10(0.01),lty=3,col="black",lwd=0.8)
  p1<-p1+ggtitle(label)+theme(plot.title = element_text(hjust = 0.5,size=30,face="bold"))+theme(axis.title =element_text(size=20))
  print(p1)
dev.off()
})
}

plot_mutideg <- function(df,incluster='cluster',infc='avg_log2FC',ingene='gene',pos=TRUE,prefix='out',
                         mycol=NULL){
df$cluster=unlist(df[incluster])
df$avg_log2FC=unlist(df[infc])
df$gene=unlist(df[ingene])
df$label <- ifelse(df$p_val_adj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")
top10sig <- df %>% group_by(cluster) %>% top_n(10,wt=avg_log2FC)
df$size <- case_when(!(df$gene %in% top10sig$gene)~ 1,
                     df$gene %in% top10sig$gene ~ 2)
dt <- filter(df,size==1)

b1 <- top10sig %>% group_by(cluster) %>% summarize(bar = max(avg_log2FC,na.rm = T))
b2 <- top10sig %>% group_by(cluster) %>% summarize(bar = min(avg_log2FC,na.rm = T))

lsed= c(1:length(unique(df$cluster)))-1
dfbar<-data.frame(x=lsed,
                  y=b1$bar)
dfbar1<-data.frame(x=lsed,
                   y=b2$bar)
iny=0
if (pos){
p1 <- ggplot()+
geom_col(data = dfbar,
         mapping = aes(x = x,y = y),
         fill = "#dcdcdc",alpha = 0.6)
iny=-0.5
}else{
p1 <- ggplot()+
geom_col(data = dfbar,
         mapping = aes(x = x,y = y),
         fill = "#dcdcdc",alpha = 0.6)+
geom_col(data = dfbar1,
         mapping = aes(x = x,y = y),
         fill = "#dcdcdc",alpha = 0.6)
}

p2 <- p1+
geom_jitter(data = dt,
            aes(x = cluster, y = avg_log2FC, color = label),
            size = 0.85,
            width =0.4)+
geom_jitter(data = data.frame(top10sig),
            aes(x = cluster, y = avg_log2FC, color = label),
            size = 1,
            width =0.4,
            shape=17)

dfcol<-data.frame(x=lsed,
                  y=iny,
                  label=unique(df$cluster))
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.4,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)

p4 <- p3+
geom_text_repel(
                data=top10sig,
                aes(x=cluster,y=avg_log2FC, label=gene),
                force = 1.2,
                arrow = arrow(length = unit(0.008, "npc"),
                              type = "open", ends = "last")
                )

p5 <- p4 +
scale_color_manual(name=NULL,
                   values = c("red","black"))

p6 <- p5+
labs(x="Cluster",y="average logFC")+
geom_text(data=dfcol,
          aes(x=x,y=y,label=label),
          size =6,
          color ="white")
p7 <- p6+
theme_minimal()+
theme(
      axis.title = element_text(size = 13,
                                color = "black",
                                face = "bold"),
      axis.line.y = element_line(color = "black",
                                 size = 1.2),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = c(1,0),
      legend.text = element_text(size = 15)
      )

pdf(paste(prefix,'_multideg.pdf'),18,10)
print(p7)
dev.off()
}

get_stat_pair <- function(inlist) {
ltime1 <- inlist
my_comparisons1 <- list()
n=1
nlen <- length(ltime1)
for (base1 in ltime1) {
c1 <- which(ltime1 %in% base1)
ltime1 <- ltime1[-c1]
for (tp in ltime1) {
my_comparisons1[[n]] <- c(base1,tp)
n <- n+1
}
}
return(my_comparisons1)
}

sdimplot_out <- function(obj,label='umap',corl=cor,pt.size=1.5,group1='sub',group2=NULL,
                         hei=10,raster=F){
obj$sub <- obj@meta.data[,group1]
rat <- cal_wvsh(obj)
wei=rat*hei
pdf(paste0(label,'_',group1,'.pdf'),wei,hei)
p1 <- DimPlot(obj, reduction = "umap", label = T, label.size = 0,repel =T,group.by=group1,
              pt.size=pt.size,raster=raster, order= sort(unique(obj$sub),decreasing=T))+scale_color_manual(values=corl)+NoAxes()
p2 <- DimPlot(obj, reduction = "umap", label = T, label.size = 0,repel =T,group.by=group1,
              pt.size=pt.size,raster=raster, order= sort(unique(obj$sub),decreasing=T))+scale_color_manual(values=corl)+NoAxes()+NoLegend()+ggtitle('')
print(p1)
print(p2)
dev.off()
png(paste0(label,'_',group1,'.png'),wei*100,hei*100)
print(p2)
dev.off()
    
if (!is.null(group2)){
obj$major <- obj@meta.data[,group2]
pdf(paste0(label,'_',group1,'.pdf'),wei,hei)
p3 <- DimPlot(obj, reduction = "umap", label = T, group.by = group2, label.size = 0,repel =T,
              pt.size=pt.size,raster=raster, order= sort(unique(obj$major),decreasing=T))+scale_color_manual(values=corl)+NoAxes()
p4 <- DimPlot(obj, reduction = "umap", label = T, group.by = group2, label.size = 0,repel =T,
              pt.size=pt.size,raster=raster, order= sort(unique(obj$major),decreasing=T))+scale_color_manual(values=corl)+NoAxes()+NoLegend()+ggtitle('')
print(p3)
print(p4)
dev.off()
png(paste0(label,'_',group2,'.png'),wei/100,hei/100)
print(p4)
dev.off()
}
}


cal_wvsh <- function(obj){
wid <- max(obj@reductions$umap@cell.embeddings[,1])-min(obj@reductions$umap@cell.embeddings[,1])
hei <- max(obj@reductions$umap@cell.embeddings[,2])-min(obj@reductions$umap@cell.embeddings[,2])
rat <- wid/hei
return(rat)
}

splot_feature <- function(obj,lgene,sub=NULL,psize=0.5,output='./',assay='Spatial',
                          hei=10,incolor=colorRampPalette(brewer.pal(9,"Blues"))(100)[-c(1:10)],raster=NULL){
all <- obj
n1 <- dim(all)[2]
if (!is.null(sub)) {
library(stringr)
sub1 <- unlist(str_split(sub,'-'))
sub <- sub1[-1]
sub1 <- sub1[1]
c1 <- which(all@meta.data[,sub1] %in% sub)
all <- subset(all,cells=Cells(all)[c1])
}
n2 <- dim(all)[2]
DefaultAssay(all) <- assay

fp_plot <- function(obj,gene.list=NULL,output="test",raster=NULL,color=blu,pt.size=1) {
all <- obj
gene <- unique(gene.list)
c <- which(gene %in% rownames(obj))
gene <- gene[c]
for (i in 1:length(gene)){
pf <- paste0(i,"_",gene[i],"_Featureplot.pdf")
pf <- gsub("/","_",pf)
rat <- cal_wvsh(all)
wei=rat*hei
pdf(paste0(output,'_',pf),wei,hei)
if (is.null(raster)) {
p1=FeaturePlot(obj,features = gene[i],pt.size=pt.size)+scale_color_gradientn(colours = color)+NoAxes()
}else{
p1=FeaturePlot(obj,features = gene[i], raster=F,pt.size=pt.size)+scale_color_gradientn(colours = color)+NoAxes()
}
if (is.null(raster)) {
p2=FeaturePlot(obj,features = gene[i],pt.size=pt.size)+scale_color_gradientn(colours = color)+NoAxes()+NoLegend()
}else{
p2=FeaturePlot(obj,features = gene[i], raster=F,pt.size=pt.size)+scale_color_gradientn(colours = color)+NoAxes()+NoLegend()
}
print(p1)
print(p2)
dev.off()
pf1=gsub('pdf','png',pf)
png(paste0(output,'_',pf1),wei*100,hei*100)
print(p2)
dev.off()
}
}
fp_plot(all,gene.list=lgene,output=output,color=incolor,pt.size=psize)
}

ch_coor <- function(ob1,u1=NULL,u2=NULL,bk=F){
if (bk){
ob1@reductions$bk_umap=ob1@reductions$umap
}
ob1@reductions$umap@cell.embeddings[,1]=u1
ob1@reductions$umap@cell.embeddings[,2]=u2
return(ob1)
}

plot_boxplot <- function(df,x,y, incolor=cor,label='out',my_comparisons=NULL,color.by=NULL,
                        bxp.errorbar.width = 0.5,width = 0.6,size=1,method="t.test",
                         p.type=c('p.signif',"p.format")[1]) {
if (is.null(my_comparisons)) {
my_comparisons=get_stat_pair(unique(df[,x]))
}
if (is.null(color.by)){
color.by=x
}
pb1 <-ggboxplot(df,x=x,y=y,color=x,bxp.errorbar.width = bxp.errorbar.width,width = width, size=size)+
      theme_classic() +
      theme(panel.grid.major=element_line(colour=NA),
      plot.title = element_text(hjust = 0.5,size=20),
      axis.title.x=element_text(size=25),
      axis.title.y=element_text(size=25),
      axis.text=element_text(size=24,face = "bold"),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
      legend.text=element_text(colour="black",size=24),
      legend.title=element_text(colour="black",size=25)) +
      labs(x = x, y =y, title = label) +
      scale_colour_manual(values = incolor)

pb1 <- pb1 + stat_compare_means(method=method,hide.ns = F,comparisons =my_comparisons,label=p.type)
pdf(paste0(label,'_boxplot_test.pdf'))
print(pb1)
dev.off()
}
get_UCell <- function(obj,inpath,ncores=15,downsample=500){
DefaultAssay(obj)='RNA'
ltf=readLines(paste0(inpath,'/01.tf/use_gene.list'))
lhvg=readLines(paste0(inpath,'/02.ntf/use_gene.list'))
shvg=readLines(paste0(inpath,'/03.ntf_sm/use_gene.list'))
markers=list('TFs'=ltf,'HVGs'=lhvg,'sHVGs'=shvg)
Idents(obj)='ingroup'
obj=subset(obj,downsample=downsample)
obj <- AddModuleScore_UCell(obj, features = markers,ncores=ncores)
return(obj)
}

get_UCellexp_plot <- function(obj1,obj2,inpath,label='out',
                         features=c('TFs_UCell','HVGs_UCell','sHVGs_UCell'),
                        ingoup='ingroup'){
inpath=inpath
plot_UCell <- function (obj, assay = NULL, features = c("nFeature_RNA", "nCount_RNA", 
    "percent.mt"), group.by = "seurat_clusters", pt.size = 0, 
    pf = NULL, wid = 1800, hei = 1000, raster = NULL, stack = TRUE, 
    flip = TRUE){
p1 <- VlnPlot(obj, features = features, group.by = group.by, 
        pt.size = pt.size, assay = assay, stack = stack, flip = flip)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
    ggtitle(pf)+NoLegend()
return(p1)
}
ob1=get_UCell(obj=obj1,inpath=inpath,ncores=15,downsample=500)
ob2=get_UCell(obj=obj2,inpath=inpath,ncores=15,downsample=500)
p1=plot_UCell(ob1,features = features,group.by = ingoup,pf = 'ob1')
p2=plot_UCell(ob2,features = features,group.by = ingoup,pf = 'ob2')
pdf(paste0(label,'_UCell.pdf'),5,5)
print(p1/p2)
dev.off()
}

plot_venn <- function(data,artist,word,incolor,label='out',size=1){
data[,'artist']=data[,artist]
data[,'word']=data[,word]
lname=as.vector(unique(data$artist))
p1=venn.diagram(
  x = list(
    data %>% filter(artist==lname[1]) %>% select(word) %>% unlist() ,
    data %>% filter(artist==lname[2]) %>% select(word) %>% unlist() ,
    data %>% filter(artist==lname[3]) %>% select(word) %>% unlist()
    ),
  category.names = lname,output = FALSE,
          filename=NULL,
          lwd = 1,
          col=incolor,
          fill = c(alpha(incolor[1],0.3), alpha(incolor[2],0.3), alpha(incolor[3],0.3)),
          cex = 0.5,
          fontfamily = "sans",
          cat.cex = 0.6,
          cat.default.pos = "outer",
          cat.pos = c(-27, 27, 135),
          cat.dist = c(0.055, 0.055, 0.085),
          cat.fontfamily = "sans",
          cat.col = incolor,
          rotation = 1
        )
pdf(paste0('venn_',label,'.pdf'),4,4)
grid.draw(p1)
dev.off()
}

