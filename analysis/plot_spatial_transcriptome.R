source("functions/function_merge.R")
source('functions/function_integ.R')
source('functions/function_plot.R')

#teleost
spe2 = Read10X('/data/work/revision/14.revision/06.validation_datasets/01.fish/B2A/out/')
image2 <- Read10X_Image('/data/work/revision/14.revision/06.validation_datasets/01.fish/B2A/spatial//sp_data//JTS12_296_A1_fr//outs//spatial')
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")
image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
spe2=seob_cluster(spe2,is.subset =F)
spe1=ch_coor(spe2,u1 = spe2@images$slice1@coordinates$row,u2=spe2@images$slice1@coordinates$col,bk = T)
spe3=subset(spe1,cells=Cells(spe1)[which(spe1@reductions$umap@cell.embeddings[,2]>50)])
spe3$ingroup='Pallium'
spe3$ingroup[which(spe3$seurat_clusters %in% c(2,8,3,7))]='SubPallium'
sdimplot_out(spe3,label='fish',corl=c('lightgreen','darkgreen'),pt.size=12,group1='ingroup',group2=NULL,hei=10)
output='./fish/'
dir.create(output)
lgene=c('thrb','nr3c1')
splot_feature(spe3,lgene,sub=NULL,psize=10,output=output,assay='Spatial',
                          hei=10,incolor=colorRampPalette(brewer.pal(9,"Blues"))(100)[-c(1:10)])

#macaque
ma=readRDS('/data/work/h1/03.figures/01.spatial/01.rds/macaca_spatial_20221101.rds')
ma$ingroup='Pallium'
ma$ingroup[which(ma$region3 %in% c('Striatum_Cd','Striatum_Pu'))]='SubPallium'
sdimplot_out(ma,label='ma',corl=c('lightgreen','darkgreen'),pt.size=0.1,group1='ingroup',group2=NULL,hei=10)
output='./ma/'
dir.create(output)
lgene=c('THRB','NR3C1','STAT6','DLX2','RARB')
splot_feature(ma,lgene,sub=NULL,psize=0.1,output=output,assay='Spatial',
                          hei=10,incolor=colorRampPalette(brewer.pal(9,"Blues"))(100)[-c(1:10)])
