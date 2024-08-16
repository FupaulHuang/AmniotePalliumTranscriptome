
library(Seurat)
library(data.table)
library(dplyr)

condense_cell <- function(obj,group.by=NULL,size=NULL,slot='data',assay='RNA',fun=mean){
mat <- data.table(t(as.matrix(slot(obj@assays[[assay]],slot))),keep.rownames=T)
meta <- obj@meta.data
meta$Celltype <- meta[,group.by]

cellist <- rownames(meta)

lcod <- c()
lscell <- c()
for (i in unique(meta$Celltype)) {
c1 <- which(meta$Celltype==i)
cellist1 <- cellist[c1]
cod <- floor(seq_along(cellist1)/size)
t1 <- data.frame(table(cod))
t1[,1] <- as.vector(t1[,1])
t1[,2] <- as.vector(t1[,2])
c2 <- which(t1[,2]<size)
short <-t1[,1][c2]
c3 <- which(cod %in% short)
cod[c3] <- t1[,1][which(t1[,2]>=size)[1]]
cod <- paste0(i,"_Cell",cod)
scell <- sample(cellist1)
lcod <- c(lcod,cod)
lscell <- c(lscell,scell)
}

map_df <- data.frame(lscell,lcod)
colnames(map_df) <- c("old_id","new_id")
rownames(map_df) <- map_df$old_id
map_df <- map_df[mat$rn,]
mat$rn <- map_df$new_id

mat <- mat[,lapply(.SD,fun), by=rn]
rownames(mat) <- mat$rn
gename <- colnames(mat)
celname <- rownames(mat)
mat <- t(mat[,rn:=NULL])
colnames(mat) <- celname
rownames(mat) <- gename
obj <- CreateSeuratObject(mat)
return(obj)
}

library(igraph)

plot_network <- function(edges,nodes=NULL,group.by=NULL,weight=NULL,coul = NULL,layout=layout.sphere,pf=NULL,main=""){
if (!is.null(weight)) {
edges$weight <- edges[,weight]
}else{
edges$weight <- 1
}
if (!is.null(group.by)) {
nodes$group <- nodes[,group.by]
}else{
nodes$group <- 'group1'
}
network <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)
library(RColorBrewer)
nlen <- length(unique(V(network)$group))

if (is.null(coul)) {
library(RColorBrewer)
coul  <- c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"))
}
coul <- coul[1:nlen]
vertex.color <- coul[as.numeric(as.factor(V(network)$group))]

pdf(paste0(pf,'_',group.by,'_network.pdf'),18,10)
plot(network, vertex.color=vertex.color, edge.width=E(network)$weight*2,layout=layout,main=main)
legend("topright", legend=levels(as.factor(V(network)$group)), 
       col = coul, bty = "n", pch=20 , pt.cex = 3, 
       cex = 1.5, text.col=coul , horiz = FALSE, 
       inset = c(0.1, 0.1))
dev.off()
}

##################################################
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

run_wgcna <- function(mat,type = "unsigned",corType = "pearson", pct.mad=0.75, label='test',trait=NULL,tom=F,RsquaredCut=0.85,hug.top=50,power=NULL){
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)

dataExpr <- mat

#Data QC
##################################################
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- as.matrix(dataExpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 1-pct.mad))[2],0.01)),])
dataExpr <- as.data.frame(t(dataExprVar))

gsg = goodSamplesGenes(dataExpr, verbose = 3)

if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

sampleTree = hclust(dist(dataExpr), method = "average")
pdf('fig1.sampleTree.pdf',18,10)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
##################################################

#calculate power
##################################################
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,networkType=type, verbose=5, RsquaredCut=RsquaredCut)

pdf('fig2.powers.pdf',18,10)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
dev.off()

saveRDS(sft,'sft.rds')
saveRDS(dataExpr,'dataExpr.rds')
if (is.null(power)) {
power = sft$powerEstimate
}
print(paste0('The power value is : ',power))
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
          ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
          ifelse(type == "unsigned", 6, 12))       
          )
          )
}
##################################################

#bulid co-expression network
##################################################
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0(label, ".tom"),
                       verbose = 3)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

pdf('fig3.plotDendroAndColors.pdf',18,10)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

try({
pdf('fig4.plotEigengeneNetworks.pdf',18,10)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()
})
##################################################
#TOMplot
#TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
plotTOM = dissTOM^7
diag(plotTOM) = NA

if (tom) {
pdf('fig5.TOMplot.pdf',18,10)
TOMplot(plotTOM, net$dendrograms, moduleColors, 
                main = "Network heatmap plot, all genes")
dev.off()
}
##################################################
#export network information
##################################################
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)

cyt = exportNetworkToCytoscape(TOM,
             edgeFile = paste(label, ".edges.txt", sep=""),
             nodeFile = paste(label, ".nodes.txt", sep=""),
             weighted = TRUE, threshold = 0,
             nodeNames = probes, nodeAttr = moduleColors)

edges <- read.table(paste(label, ".edges.txt", sep=""),sep='\t',header=T)
nodes <- read.table(paste(label, ".nodes.txt", sep=""),sep='\t',header=T)

library(RColorBrewer)
coul<-c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"))

colnames(nodes) <- c(colnames(nodes)[1:2],'group')
lgroup <- unique(nodes$group)
lgroup <- lgroup[-grep('grey',lgroup)]
for (i in lgroup) {
nodes1 <- subset(nodes,group==i)
c1 <- which(edges[,1] %in% nodes1[,1])
edges1 <- edges[c1,]
c1 <- which(edges1[,2] %in% nodes1[,1])
edges1 <- edges1[c1,]
mtitle <- dim(nodes1)[1]
plot_network(edges=edges1,nodes=nodes1,group.by=NULL,weight='weight',coul = coul,layout=layout.sphere,pf=paste0("all_genes_",i),main=mtitle)
}

##################################################
#get hub genes
##################################################

library(dplyr)
edgeData1 <- cyt$edgeData[,c(1,2,3)]
edgeData2 <- cyt$edgeData[,c(2,1,3)]
nodeattrib <- cyt$nodeData[,c(1,3)]
colnames(nodeattrib) <- c("nodeName", "Module")
colnames(edgeData1) <- c("Node1","Node2","Weight")
colnames(edgeData2) <- c("Node1","Node2","Weight")
edgeData <- rbind(edgeData1, edgeData2)
edgeData$Module1 <- nodeattrib[match(edgeData$Node1, nodeattrib$nodeName), 2]
edgeData$Module2 <- nodeattrib[match(edgeData$Node2, nodeattrib$nodeName), 2]
edgeData <- edgeData[edgeData$Module1==edgeData$Module2,c(1,3,4)]

nodeTotalWeight <- edgeData %>% group_by(Node1, Module1) %>% summarise(weight=sum(Weight))
nodeTotalWeight <- nodeTotalWeight[with(nodeTotalWeight, order(Module1, -weight)),]
nodeTotalWeightTop = nodeTotalWeight %>% group_by(Module1) %>% top_n(hug.top, weight) %>% data.frame()
write.table(nodeTotalWeightTop,'hug_gene_list.xls',sep='\t',quote=F)
lgroup <- unique(nodeTotalWeightTop$Module1)
lgroup <- lgroup[-grep('grey',lgroup)]

for (i in lgroup) {
nodes1 <- subset(nodeTotalWeightTop,Module1==i)
c1 <- which(edges[,1] %in% nodes1[,1])
edges1 <- edges[c1,]
c1 <- which(edges1[,2] %in% nodes1[,1])
edges1 <- edges1[c1,]

c1 <- which(nodes[,1] %in% nodes1[,1])
nodes1 <- nodes[c1,]
nodes1 <- subset(nodes1,group==i)

mtitle <- dim(nodes1)[1]
plot_network(edges=edges1,nodes=nodes1,group.by=NULL,weight='weight',coul = coul,layout=layout.sphere,pf=paste0("hug_genes_",i),main=mtitle)
}


#calculate correlation between gene modules and trait

if(!is.null(trait)) {
  traitData = trait
  sampleName = rownames(dataExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]

if (corType=="pearsoon") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf('fig6.labeledHeatmap.pdf',18,10)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
}
##################################################
#save result
ob1 <- list(corType,dataExpr, MEs_col, traitData, nSamples, net)
save(corType,dataExpr, MEs_col, traitData, nSamples, net, file=paste0(label,'_result.RData'))
return(ob1)
}

cor_mt <- function(indata,module,pheno){
load(indata)
if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
             as.matrix(geneModuleMembership), nSamples))
}else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
             as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
moduleGenes = moduleColors == module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))

pdf(paste0(module,"_",pheno,'_verboseScatterplot.pdf'),18,10)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
}


plot_heatmap <- function(trait,mergedMEs,pf='test'){
library(pheatmap)
cor_ADR <- signif(WGCNA::cor(trait,mergedMEs,use="p",method="pearson"),5)
p.values <- corPvalueStudent(cor_ADR,nSamples=nrow(trait))
pdf(paste0(pf,'_trait_Heatmap.pdf'),18,10)
pheatmap(cor_ADR,display_numbers = matrix(ifelse(p.values <= 0.01, "**", ifelse(p.values<= 0.05 ,"*"," ")), nrow(p.values)),fontsize=18)
dev.off()
}

get_module_gene <- function(dataExpr,module,moduleColors){
probes = colnames(dataExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule]
return(modProbes)
}

seurat_wgcna <- function(obj,group.by=NULL,size=NULL,slot='data',assay='RNA',fun=mean,
                         type = "unsigned",corType = "pearson", pct.mad=0.75, label='test',trait=NULL,nThreads=10){
enableWGCNAThreads(nThreads=nThreads)
if (!is.null(size)) {
obj <- condense_cell(obj,group.by=group.by,size=size,slot=slot,assay=assay,fun=fun)
saveRDS(obj,'condense.rds')
}
trait <- obj@meta.data
run_wgcna(mat=obj@assays$RNA@counts,type = type, corType = corType, pct.mad=pct.mad, label=label,trait=trait)
}

get_module <- function(inpath,label,ingene){
df1=read.table(paste0('01.wgcna/',inpath,'/hug_gene_list.xls'),sep='\t',header=T)
c1 <- which(df1$Node1 %in% ingene)
df2=df1[c1,]
print(df2)
c1 <- which(df1$Module1 %in% df2$Module1)
df1=df1[c1,]
df1$Module1=paste0(label,'_',df1$Module1)
return(df1)
}

get_score <- function(all,inlist,label){
DefaultAssay(all) <- 'RNA'
out <- all@meta.data[,c('first','second','third','type')]
for (i in unique(inlist$Module1)) {
df1 <- subset(inlist,Module1==i)
ingene <- df1$Node1
tmp <- AddModuleScore(all,features=list(ingene))
out[,i] <- 'unknown'
out[,i] <- tmp$Cluster1
#tmp <- AddModuleScore(all,features=list(ingene),name=xnam)
}
write.table(out,paste0(label,'_addscore.xls'),sep='\t',quote=F)
}

plot_module <- function(df1,group_by,label='out'){
df1$cluster=df1[,group_by]
df4=df1
#df4 <- df1[order(df1$cluster),]
df4 <- df4[,-c(1:4)]


ldf=list()
for(i in colnames(df4)[-dim(df4)[2]]) {
df4 %>% group_by(cluster)%>% summarise(y=mean(get(i))) -> df5
colnames(df5) <- c('cluster',i)
ldf[[i]]=df5
}

df5=Reduce(cbind,ldf)
rownames(df5)=df5$cluster
df5=df5[,-grep('cluster',colnames(df5))]
p1 <- pheatmap(df5,cluster_cols = T,cluster_rows = T,,border='black',scale="none",color =  colorRampPalette(c("white", "#7B2411"))(100),,border_color='black',cellwidth=15,cellheight=15)

pdf(paste0(label,'_',group_by,'.pdf'),6,6)
print(p1)
dev.off()
}



