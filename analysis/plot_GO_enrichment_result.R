source("functions/function_enrichment.R")

#GO enrichment using one gene set
df1<-read.table('conserved_celltype.xls',sep='\t',header=T)
ingene=rownames(df1)
test = bitr(ingene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb='org.Hs.eg.db')
ingene=as.vector(test$ENTREZID)
ob1 <- enrich_result(ingene,fun='GO',label='vip_lge_conserved')

#GO enrichment using two or more gene sets
df1 <- read.table('/hwfssz1/ST_SUPERCELLS/P22Z10200N0664/huangbaoqian/02.macaca/03.figures/02.deg/02.cluster_deg/in_deg_markers.xls',sep='\t',header=T)
c1 <- which(df1$cluster %in% c(7,8,9))
df1 <- df1[c1,]
df1 <- df1[df1$FC > 2,]
df1 <- df1[df1$p_val_adj < 0.05,]
lgene=list()
for (i in unique(df1$cluster)) {
ingene=as.vector(df1$gene[df1$cluster==i])
test = bitr(ingene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb='org.Hs.eg.db')
ingene=as.vector(test$ENTREZID)
lgene[[as.character(i)]] <- ingene
}

ob1 <- compare_result(lgene,fun='GO',label='3IN')
ob1=filter(ob1, p.adjust < 0.05)
c1=grep('epithel',ob1@compareClusterResult$Description)
ob2=filter(ob1, ob1@compareClusterResult$ID %in% ob1@compareClusterResult$ID[-c1])
lgo <-c(
'GO:0030900',
'GO:0010763',
'GO:0021537',
'GO:0061351',
'GO:0051386',
'GO:0060828',
'GO:0038007',
'GO:0007405',
'GO:0048015',
'GO:0048017')
ob2=filter(ob2, ob2@compareClusterResult$ID %in% lgo)
ob3=filter(ob2, p.adjust < .05, qvalue < 0.2)

compare_plot(eobj=ob2,label='3IN_subset',colours = c('#336699','#66CC66','#FFCC33'),
                       fun= "GO", showCategory = 10,
                       categorySize="pvalue",foldChange=NULL,node_label="all",color_category='firebrick',color_gene='steelblue',
                       legend_n=2,inpie="count", cex_category=1.5, layout="kk",wid=10,hei=10)

