source("functions/function_plot.R")

#Plot Major Cell Types of six species
meta <- read.table('all_metadata_20230714.xls',sep='\t',header=T)
lsub=c('IN','EX','AST','OLI', 'OPC','MIC','PERI','ENDO','VLMC','MUR','NPC','EPEN','EG')
ltype=c('IN','EX','AST','OLI', 'OPC','MIC','PERI','ENDO','VLMC','MUR','NPC','EPEN','EG')
lcor <- cor[1:length(ltype)]
names(lcor) <- ltype

#lizard
plot_umap(inrds='01.rds/03.human_homology_rds/li_human_homology.rds',
          inmeta='01.rds/lizard_metadata_20230713.xls',
          label='lizard')

#human
plot_umap(inrds='01.rds/03.human_homology_rds/hs.rds',
          inmeta='01.rds/human_metadata_20230713.xls',
          label='human')

#macaque
plot_umap(inrds='01.rds/03.human_homology_rds/ma_human_homology.rds',
          inmeta='01.rds/macaque_metadata_20230713.xls',
          label='macaque')

#mouse
plot_umap(inrds='01.rds/03.human_homology_rds/mm_human_homology.rds',
          inmeta='01.rds/mouse_metadata_20230713.xls',
          label='mouse')

#turtle
plot_umap(inrds='01.rds/03.human_homology_rds/tt_human_homology.rds',
          inmeta='01.rds/turtle_metadata_20230713.xls',
          label='turtle')

#zebrafinch
plot_umap(inrds='01.rds/03.human_homology_rds/zf_human_homology.rds',
          inmeta='01.rds/zebrafinch_metadata_20230713.xls',
          label='zebrafinch')

meta$Major <- dplyr::recode(meta$Major,'noneuron'='Non-Neuronal')
meta$Major <- dplyr::recode(meta$Major,'NPC'='Non-Neuronal')

#Bar plot of cell proportions
pct_barplot(meta,group.by='Major',split.by='order',label='major_celltype')


#Plot macaque cell subtypes
source("functions/function_merge.R")
source('functions/function_integ.R')

dimplot_out <- function(obj,label='umap',corl=cor,pt.size=1.5,group1='sub',group2='major',wei=10,hei=10){
obj$sub <- as.vector(obj@meta.data[,group1])
obj$major <- as.vector(obj@meta.data[,group2])
pdf(paste0(label,'_',group1,'_',group2,'.pdf'),wei,hei)
p1 <- DimPlot(obj, reduction = "umap", label = T, label.size = 0,repel =T,group.by=group1,pt.size=0.5,raster=F)+scale_color_manual(values=corl)+NoAxes()
p2 <- DimPlot(obj, reduction = "umap", label = T, label.size = 0,repel =T,group.by=group1,pt.size=pt.size,raster=F)+scale_color_manual(values=corl)+NoAxes()+NoLegend()+ggtitle('')
p3 <- DimPlot(obj, reduction = "umap", label = T, group.by = group2, label.size = 0,repel =T,pt.size=0.5,raster=F)+scale_color_manual(values=corl)+NoAxes()
p4 <- DimPlot(obj, reduction = "umap", label = T, group.by = group2, label.size = 0,repel =T,pt.size=pt.size,raster=F)+scale_color_manual(values=corl)+NoAxes()+NoLegend()+ggtitle('')
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()
}


cal_wvsh <- function(obj){
wid <- max(obj@reductions$umap@cell.embeddings[,1])-min(obj@reductions$umap@cell.embeddings[,1])
hei <- max(obj@reductions$umap@cell.embeddings[,2])-min(obj@reductions$umap@cell.embeddings[,2])
rat <- wid/hei
return(rat)
}

all <- readRDS('01.rds/03.human_homology_rds/ma_human_homology.rds')
all$Celltype <- all$sub
rat <- cal_wvsh(all)
cor <- cor[-5]
dimplot_out(all,label='fig1b',group1='predict_Celltype',group2='Celltype',wei=10*rat,hei=10,corl=cor,pt.size=1)

gene=c('SLC17A7','GAD1')
pmarkers(all,gene)

