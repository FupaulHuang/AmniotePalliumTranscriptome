
library(future)
plan("multicore", workers = 4)
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
library(DOSE)
library(enrichplot)


sub_deg <- function(obj, label='out', group.by=NULL,only.pos = F, return.thresh = 0.01, logfc.threshold = 0.25, min.pct = 0.1, wt='avg_log2FC'){
all <- obj
Idents(all) <- group.by
all_markers <- FindAllMarkers(all, only.pos = only.pos, return.thresh = return.thresh, logfc.threshold = logfc.threshold, verbose = F)

all_markers$FC <- 2^(all_markers[,wt])
write.table(all_markers,paste0(label,'_',group.by,'_deg_markers.xls'),sep='\t',quote=F)
return(all_markers)
}

deg_dotplot <- function(obj, all_markers, group.by=NULL, pattern=NULL,ntop=5,wt='avg_log2FC',inlevels=NULL,label='out'){
all <- obj
Idents(all) <- group.by
if (!is.null(pattern)) {
c1 <- grep(pattern,all_markers$gene)
if (length(c1)!=0) {
all_markers <- all_markers[-c1,]
}
}
#all_markers <- all_markers[order(all_markers$cluster),]
all_markers$avg_log2FC <- all_markers[,wt]
all_markers$FC <- 2^(all_markers$avg_log2FC)

top5 <- all_markers %>% group_by(cluster) %>% top_n(n=ntop, wt=avg_log2FC)

blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)

if (is.null(inlevels)) {
inlevels <- sort(levels(all))
}

all@active.ident <- factor(all@active.ident,levels=inlevels)

hei <- ceiling(length(unique(top5$gene))*0.2)
hei <- max(hei,10)
pdf(paste0(label,"_",group.by,'_deg_dotplot.pdf'),10,hei)
print(DotPlot(all,features=unique(top5$gene))+coord_flip()+scale_color_gradientn(colours = blu)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)))
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

enrich_plot <- function(eobj,label='out',colours = c('#336699','#66CC66','#FFCC33'),
                       fun= "GO", showCategory = 10,
                       categorySize="pvalue",foldChange=NULL,node_label="all",color_category='firebrick',color_gene='steelblue',
                       cex_category=1.5,layout="kk",wid=18,hei=10) {
pdf(paste0(label,"_enrich_",fun,"_plot.pdf"),wid,hei)

ttl <- paste0(label,"_",fun)
gttl <-ggtitle(ttl)  

p1 <- dotplot(eobj,showCategory = showCategory,title=ttl)+ scale_color_gradientn(values = seq(0,1,0.2),colours = colours)
p2 <- barplot(eobj, showCategory=showCategory, title=ttl)+ scale_fill_gradientn(values = seq(0,1,0.2),colours = colours)
p3 <- mutate(eobj, qscore = -log(p.adjust, base=10)) %>% barplot(x="qscore",showCategory=showCategory, title=ttl)+ scale_fill_gradientn(values = seq(0,1,0.2),colours = colours)
p4 <- cnetplot(eobj, categorySize=categorySize, foldChange=foldChange,node_label=node_label,color_category=color_category,color_gene=color_gene)+ gttl
p5 <- cnetplot(eobj, foldChange=foldChange, circular = TRUE, colorEdge = TRUE,node_label=node_label, color_category=color_category,color_gene=color_gene)+ gttl
p6 <- heatplot(eobj, foldChange=foldChange, showCategory=showCategory) + scale_color_gradientn(values = seq(0,1,0.2),colours = colours)+gttl

eobj1 <- pairwise_termsim(eobj) 

p9 <- emapplot(eobj1, cex_category=cex_category,layout=layout) + gttl+ scale_color_gradientn(values = seq(0,1,0.2),colours = colours)
p10 <- upsetplot(eobj)+ gttl


print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
try(print(p9))
print(p10)

try({
if (exists('treeplot')) {
p7 <- treeplot(eobj1)
p8 <- treeplot(eobj1, hclust_method = "average")
print(p7)
print(p8)
}
})
dev.off()
}



enrich_result <- function(vgene,p.val=0.05,OrgDb='org.Hs.eg.db',label='out',
                          keyType='ENTREZID',colours = c('#336699','#66CC66','#FFCC33'), pAdjustMethod='BH',
                       fun= "GO", q.val=0.2, ont = "BP", showCategory = 10,organism = "hsa",use_internal_data=T,
                       minGSSize= 5,maxGSSize= 500,
                       categorySize="pvalue",foldChange=NULL,node_label="all",color_category='firebrick',color_gene='steelblue',interm=NULL,
                       wid=18,hei=10){

fun.use=paste0('enrich',fun)

if (fun=="GO"){
eobj <- enrichGO(gene         = vgene,
                OrgDb         = OrgDb,
                keyType       = keyType,
                ont           = ont,
                pAdjustMethod = pAdjustMethod,
                pvalueCutoff  = p.val,
                qvalueCutoff  = q.val,
                readable=TRUE)
p11 <- goplot(eobj)
png(paste0(label,"_GO_goplot.png"),1800,1000)
print(p11)
dev.off()
}

if (fun=="KEGG"){
kk <- enrichKEGG(gene         = vgene,
                 organism     = organism,
                 pvalueCutoff = p.val,
                 use_internal_data=use_internal_data)

eobj <- setReadable(kk,OrgDb=OrgDb,keyType=keyType)

}

if (fun=="DO"){
eobj <- enrichDO(gene       = vgene,
              ont           = "DO",
              pvalueCutoff  = p.val,
              pAdjustMethod = pAdjustMethod,
              minGSSize     = minGSSize,
              maxGSSize     = maxGSSize,
              qvalueCutoff  = q.val,
              readable      = TRUE)

}
 
if (fun=="NCG"){
eobj <- enrichNCG(gene      = vgene,
              pvalueCutoff  = p.val,
              pAdjustMethod = pAdjustMethod,
              minGSSize     = minGSSize,
              maxGSSize     = maxGSSize,
              qvalueCutoff  = q.val,
              readable      = TRUE)


}

if (fun=="DGN"){
eobj <- enrichDGN(gene      = vgene,
              pvalueCutoff  = p.val,
              pAdjustMethod = pAdjustMethod,
              minGSSize     = minGSSize,
              maxGSSize     = maxGSSize,
              qvalueCutoff  = q.val,
              readable      = TRUE)

}

if (fun=="enricher"){

x <- enricher(gene      = vgene,
              pvalueCutoff  = p.val,
              pAdjustMethod = pAdjustMethod,
              minGSSize     = minGSSize,
              maxGSSize     = maxGSSize,
              qvalueCutoff  = q.val,
              TERM2GENE = interm)

eobj <- setReadable(x,OrgDb=OrgDb,keyType=keyType)
}

saveRDS(eobj, paste0(label,"_",fun,"_enrich.rds"))
out=eobj@result
write.table(out,paste0(label,"_enrich_",fun,"List.xls"),row.names = FALSE,quote = FALSE,sep = "\t")
enrich_plot(eobj=eobj,label=label,colours = colours,
                       fun= fun, showCategory = showCategory,
                       categorySize=categorySize,foldChange=foldChange,node_label=node_label,color_category=color_category,color_gene=color_gene,
                       wid=wid,hei=hei)
return(eobj)
}


filter_genes <- function(indeg,fc.filter=1.5,p.val=0.05,pct=0.1, OrgDb='org.Hs.eg.db') {
filter_gene <- indeg
filter_gene<-filter_gene[filter_gene$p_val_adj<= p.val & filter_gene$pct.1 >= pct,]
filter_gene<-filter_gene[filter_gene$FC >= fc.filter | filter_gene$FC <= 1/fc.filter,]

#df1 <- filter_gene %>% group_by(cluster) %>% sample_n(20)

c=unique(filter_gene$cluster)
n=length(unique(filter_gene$cluster))
gene=list()
symbol=list()
j =1
for(x in c){
      sub=subset(filter_gene,filter_gene$cluster==x)
  test = bitr(sub$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=OrgDb)
    gene[[j]]=as.character(test$ENTREZID)
    symbol[[j]] <- as.character(test$SYMBOL)
      j=j+1
}


names(gene)=c
names(symbol)=c
lgene <- gene
return(list(gene,symbol))
}



bardot_plot <- function(comp_result,showCategory=40,title='test',wid=18,hei=10,facet=F){
pdf(paste0(title,"_bardotplot.pdf"),wid,hei)
barplot(comp_result,showCategory = showCategory,title = title)
dotplot(comp_result,showCategory = showCategory,title = title)
if (facet){
fac <- facet_grid(ONTOLOGY~., scale="free")
barplot(comp_result,showCategory = showCategory,title = title)+fac
dotplot(comp_result,showCategory = showCategory,title = title)+fac
}

dev.off()
}
#ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

compare_plot <- function(eobj,label='out',colours = c('#336699','#66CC66','#FFCC33'),
                       fun= "GO", showCategory = 10,
                       categorySize="pvalue",foldChange=NULL,node_label="all",color_category='firebrick',color_gene='steelblue',
                       legend_n=2,inpie="count", cex_category=1.5, layout="kk",wid=18,hei=10) {
pdf(paste0(label,"_comparecluster_",fun,"_plot.pdf"),wid,hei)

ttl <- paste0(label,"_",fun)
gttl <-ggtitle(ttl)  

p1 <- dotplot(eobj,showCategory = showCategory,title=ttl)+ scale_color_gradientn(values = seq(0,1,0.2),colours = colours)
p4 <- cnetplot(eobj, categorySize=categorySize, foldChange=foldChange,node_label=node_label,color_category=color_category,color_gene=color_gene)+ gttl
p5 <- cnetplot(eobj, foldChange=foldChange, circular = TRUE, colorEdge = TRUE,node_label=node_label, color_category=color_category,color_gene=color_gene)+ gttl

eobj1 <- pairwise_termsim(eobj) 
p9 <- emapplot(eobj1, cex_category=cex_category,legend_n=legend_n,pie=inpie, layout=layout) + gttl+ scale_color_gradientn(values = seq(0,1,0.2),colours = colours)


print(p1)
print(p4)
print(p5)
try(print(p9))

dev.off()
}

compare_result <- function(lgene,p.val=0.05,OrgDb='org.Hs.eg.db',label='out',
                          keyType='ENTREZID',colours = c('#336699','#66CC66','#FFCC33'), pAdjustMethod='BH',
                       fun= "GO", q.val=0.2, ont = "BP", showCategory = 10,organism = "hsa",use_internal_data=T,
                       categorySize="pvalue",foldChange=NULL,node_label="all",color_category='firebrick',color_gene='steelblue',
                       minGSSize= 5,maxGSSize= 500,interm=NULL,wid=18,hei=10){

fun.use=paste0('enrich',fun)

if (fun=="GO"){
eobj <- compareCluster(geneCluster   = lgene,
                            fun           = fun.use,
                            pvalueCutoff  = p.val,
                            pAdjustMethod = pAdjustMethod,
                            OrgDb = OrgDb,
                            ont = ont,
                            readable = TRUE)
}

if (fun=="KEGG"){
kk <- compareCluster(geneCluster = lgene, 
                     fun = fun.use,
                     pvalueCutoff  = p.val,
                     pAdjustMethod = pAdjustMethod,
                     organism     = organism,
                     use_internal_data=use_internal_data
                     )

eobj <- setReadable(kk,OrgDb=OrgDb,keyType=keyType)

}

if (fun=="DO"){

eobj <- compareCluster(geneCluster = lgene, 
              ont           = "DO",
              fun = fun.use,
              pvalueCutoff  = p.val,
              pAdjustMethod = pAdjustMethod,
              minGSSize     = minGSSize,
              maxGSSize     = maxGSSize,
              qvalueCutoff  = q.val,
              readable      = TRUE)

}
 
if (fun=="NCG"){
eobj <- compareCluster(geneCluster = lgene, 
              fun = fun.use,
              pvalueCutoff  = p.val,
              pAdjustMethod = pAdjustMethod,
              minGSSize     = minGSSize,
              maxGSSize     = maxGSSize,
              qvalueCutoff  = q.val,
              readable      = TRUE)


}

if (fun=="DGN"){
eobj <- compareCluster(geneCluster = lgene, 
              fun = fun.use,
              pvalueCutoff  = p.val,
              pAdjustMethod = pAdjustMethod,
              minGSSize     = minGSSize,
              maxGSSize     = maxGSSize,
              qvalueCutoff  = q.val,
              readable      = TRUE)

}

if (fun=="enricher"){
x <- compareCluster(geneCluster = lgene, 
              fun = fun,
              pvalueCutoff  = p.val,
              pAdjustMethod = pAdjustMethod,
              minGSSize     = minGSSize,
              maxGSSize     = maxGSSize,
              qvalueCutoff  = q.val,
              TERM2GENE = interm)

eobj <- setReadable(x,OrgDb=OrgDb,keyType=keyType)
}

saveRDS(eobj, paste0(label,"_comparecluster_",fun,".rds"))

out=eobj@compareClusterResult
write.table(out,paste0(label,"_comparecluster_",fun,"List.xls"),row.names = FALSE,quote = FALSE,sep = "\t")

compare_plot(eobj=eobj,label=label,colours = colours,
                       fun= fun, showCategory = showCategory,
                       categorySize=categorySize,foldChange=foldChange,node_label=node_label,color_category=color_category,color_gene=color_gene,
                       wid=wid,hei=hei)
return(eobj)
}



