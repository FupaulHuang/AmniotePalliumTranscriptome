
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(RColorBrewer)
library(ggpubr)
library(ggsci)
library(pheatmap)

colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))
len <- 100
col<-c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"),colpalettes,rainbow(len))
cor <- col



obj_heatmap <- function(obj,group1,group2,lgene,pf='test.jpg',
                        inassay='RNA',insep=' ',ann_colors=NULL,anno_col=NULL,incolor=c('lightgrey','black'),
                        widsize=NULL,heisize=NULL,lorder=NULL,inave=NULL,scale="none"){
obj$ingroup <- paste0(obj@meta.data[,group1],insep,obj@meta.data[,group2])
obj$ingroup <- factor(obj$ingroup,levels=unique(obj$ingroup))
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
if (is.null(inave)) {
ave <- AverageExpression(obj,assays=inassay,group.by='ingroup',features=lgene)
saveRDS(ave,'ave.rds')
}
heat_df <- ave[[inassay]]
if (!is.null(lorder)) {
heat_df <- heat_df[,lorder]
}

if (is.null(anno_col)) {
anno_col <- data.frame(colnames(heat_df))
anno_col[,group1] <- gsub(paste0(insep,'.*'),'',anno_col[,1])
anno_col[,group2] <- gsub(paste0('^.*',insep),'',anno_col[,1])

rownames(anno_col) <-  colnames(heat_df)
anno_col <- anno_col[,-1]
anno_col <- arrange(anno_col,group1,group2)
}

if(is.null(ann_colors)) {
col2=col
col1=pal_npg("nrc")(10) 
col1 <- col1[1:length(unique(obj@meta.data[,group1]))] 
names(col1) <-unique(obj@meta.data[,group1])  
col2 <- col2[1:length(unique(obj@meta.data[,group2]))]
names(col2) <- unique(obj@meta.data[,group2])

ann_colors<- list(group1 = col1, group2= col2)
names(ann_colors) <- c(group1,group2)
}



if (is.null(widsize)){
widsize= max(0.3*dim(heat_df)[2]+2,10)
heisize= max(0.3*dim(heat_df)[1]+2,10)
}
heat_df=scale(heat_df)
saveRDS(heat_df,'data_heatmap.rds')
pheatmap(heat_df,
         cluster_cols = F, cluster_rows = F,
         scale=scale,color =  colorRampPalette(incolor)(100),
         filename=pf,width = widsize, height = heisize,
         annotation_col = anno_col,show_colnames=F,annotation_colors = ann_colors,annotation_names_col=T,
         cellheight=15,cellwidth=15,
         border=FALSE)
}


#obj_heatmap(obj=obj,group1='position2',group2='cluster_int_sub2',lgene=rownames(obj)[1:20],pf='test.jpg',
#                        inassay='RNA',insep=' ',ann_colors=NULL,anno_col=NULL,incolor=c('lightgrey','black'),
#                        widsize=NULL,heisize=NULL)
