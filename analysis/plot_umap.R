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

#Plot marker FeaturesPlot
gene=c('SLC17A7','GAD1')
pmarkers(all,gene)

#Plot heatmap
source('functions/function_heatmap.R')
all=readRDS('all_anno_20240717.rds')
ave=readRDS('ave.rds')
lsub=c('IN','EX','AST','OLI', 'OPC','MIC','PERI','ENDO','VLMC','MUR','NPC','EPEN','EG')
lgene <- c('GAD1','GAD2','SLC17A7','SLC17A6', 'SLIT3','SLC1A2','SLC1A3','GFAP',
           #'MOG','OPALIN','OLIG1', 
           'MBP','PLP1','BCAS1','SOX6',
           'PDGFRA',
           'PCDH15','APBB1IP',
           #'CTSS',
           'PTPRC','C1QB',
           'C1QC','PDGFRB','KCNJ8',
           'FLT1','CLDN5',
           'DCN','COL1A2','COL1A1',
           'LUM','SOX11','SOX4',
           #'SPEF2',
           'MKI67','TOP2A',
           'FOXJ1',
           'SOX9')
meta$Sub <- factor(meta$Sub,levels=lsub)
meta <- dplyr::arrange(meta,Sub,order)
meta$ingroup <- paste0(meta$Sub,' ',meta$order)
lorder <- unique(meta$ingroup)
meta=all@meta.data
group1 <- 'Sub'
group2 <- 'order'
obj <- all
col1= col
col2= pal_simpsons("springfield")(16)
col1 <- col1[1:length(unique(obj@meta.data[,group1]))]
names(col1) <-unique(obj@meta.data[,group1])
col2 <- col2[1:length(unique(obj@meta.data[,group2]))]
names(col2) <- unique(obj@meta.data[,group2])
ann_colors<- list(group1 = col1, group2= col2)
names(ann_colors) <- c(group1,group2)
obj_heatmap(obj=all,group1='Sub',group2='order',lgene=lgene,pf='markers_heatmap.pdf',
                        inassay='RNA',insep=' ',ann_colors=ann_colors,anno_col=NULL,incolor=c('white','lightgrey','black'),
                        widsize=NULL,heisize=NULL,lorder=lorder,scale='row',inave=NULL)

#Plot markers volin plot
plot_vlnplot <- function(ob1,gene,ingroup='seurat_clusters',inlabel='out',inassay='RNA',slot = "data",corl=c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"),colpalettes,rainbow(20))){
    DefaultAssay(ob1) <- 'RNA'
    gene=unique(gene)
    #blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)
    p1=VlnPlot(ob1,features = gene,group.by=ingroup,stack = T,flip =T,fill.by ='ident',slot = slot)+
            scale_fill_manual(values=corl)+
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
            NoLegend()+
            ggtitle(inlabel)

    inwid=length(unique(ob1@meta.data[,ingroup]))*28/100
    inlen=max(c(19*length(gene),884))/100
    pdf(paste0(inlabel,'_markers_vlnplot.pdf'),inwid,inlen)
    print(p1)
    dev.off()
}

markers <- read.csv(header = T, text = '
celltype,marker
L1/2_IT,SLC17A7
L1/2_IT,GPR83
L1/2_IT,PDZD2
L1/2_IT,NTNG1
L3/4_IT,MYLK
L3/4_IT,PLCH1
L3/4_IT,RORB
L3/4_IT,IL1RAPL2
L3/4_IT,OPRK1
L3/4_IT,SMYD1
L3/4_IT,ETV1
L5/6 IT Car3,TLE4
L5/6 NP,NXPH2
L5 ET,FEZF2
L5 ET,ABO
L6b,HCRTR2
L6b,SEMA3E
L6b,AMOTL1
LAMP5,GAD1
LAMP5,GAD2
LAMP5,RELN
LAMP5,ADARB2
LAMP5,EYA4
LAMP5,GAD1
LAMP5,LAMP5
PAX6,PAX6
PVALB,PVALB
SST,SST
VIP,CALB2
VIP,VIP
Astrocyte,GFAP
Astrocyte,SLC1A2
Astrocyte,SLC1A3
Oligodendrocyte,MOG
Oligodendrocyte,PLP1
OPC,PDGFRA
OPC,COL9A1
Microglia,APBB1IP
Microglia,PTPRC
Endothelial,FLT1
Endothelial,RGS5
Endothelial,COL1A2
Endothelial,PDGFRB
')
all=readRDS('/data/work/graduation/09.summary/06.pfc/macaca_pfc_20230215_diet.rds')
plot_vlnplot(all,markers$marker,ingroup='Celltype',inlabel='out',inassay='RNA',slot = "data")