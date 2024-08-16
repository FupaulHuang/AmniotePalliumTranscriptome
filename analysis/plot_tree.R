#Plot tree using HVGs
source('functions/function_integ.R')
source("functions/function_merge.R")
source("functions/function_plot.R")

library(ggtree)
library(ggplot2)
all <- readRDS('01.rds/all.rds')
meta <- read.table('01.rds/all_metadata_20230821.xls',sep='\t',header=T)
all@meta.data <- meta
all$ingroup=paste0(all$order,'_',all$Sub)

ave=AverageExpression(all,group.by='ingroup',assays ='RNA',slot ='data')

get_hvgs <- function(obj,nfeatures=2000){
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
hvg=VariableFeatures(obj)
write.table(hvg,paste0('hvg_',nfeatures,'.txt'),sep='\t',quote=F,row.names=F,col.names=F)
}
get_hvgs(all,nfeatures=3000)

lsub=c('IN','EX','AST','OLI', 'OPC','MIC','PERI','ENDO','VLMC','MUR','NPC','EPEN','EG')
ltype=c('IN','EX','AST','OLI', 'OPC','MIC','PERI','ENDO','VLMC','MUR','NPC','EPEN','EG')
lcor <- cor[1:length(ltype)]
names(lcor) <- ltype

lct=c('ENDO','MUR','PERI','VLMC')
inpattern=paste0(lct,collapse='|')
c1=grep(inpattern,colnames(ave))
inct=as.vector(colnames(ave)[-c1])

anno=data.frame(colnames(ave))
colnames(anno)='taxa'
rownames(anno)=anno[,1]
inlist=strsplit(anno[,1],'_')
anno$species=as.vector(sapply(inlist,'[',1))
anno$celltype=as.vector(sapply(inlist,'[',2))
source("/hwfssz5/ST_SUPERCELLS/P22Z10200N0639/huangbaoqian/03.graduate_20230712/02.F1/function_plot.R")
lsub=c('IN','EX','AST','OLI', 'OPC','MIC','PERI','ENDO','VLMC','MUR','NPC','EPEN','EG')
ltype=c('IN','EX','AST','OLI', 'OPC','MIC','PERI','ENDO','VLMC','MUR','NPC','EPEN','EG')
lcor <- cor[1:length(ltype)]
names(lcor) <- ltype

d1=dist(x=t(x=ave[hvg,inct]))
t1=hclust(d1)
t2=ape::as.phylo(t1)

dd=anno[inct,]
p1=ggtree(t2,branch.length='none',layout = "circular",size=2)
p1 <- p1 %<+% dd + geom_tiplab(aes(color=celltype),size=12,offset = 2)+
       geom_tippoint(aes(color=celltype),size=8)

incolor=lcor[unique(dd$celltype)]
p1=p1+scale_colour_manual(values=incolor)
p1=p1+xlim(NA, 40)
p1=p1+theme(legend.position="right")
p2=p1

p1=ggtree(t2,branch.length='none',size=2)
p1 <- p1 %<+% dd + geom_tiplab(aes(color=celltype),size=12,offset = 1)+
       geom_tippoint(aes(color=celltype),size=8)
#p1 <- p1 %<+% dd + geom_text(aes(color=place, label=label), hjust=-0.5) +
#       geom_tippoint(aes(size=value, shape=place, color=place), alpha=0.25)
incolor=lcor[unique(dd$celltype)]
p1=p1+scale_colour_manual(values=incolor)
p1=p1+xlim(NA, 70)
p1=p1+theme(legend.position="right")
pdf('tree.pdf',18,18)
print(p1)
print(p2)
dev.off()
