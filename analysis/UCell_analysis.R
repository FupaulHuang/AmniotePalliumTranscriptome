source("functions/function_merge.R")
source('functions/function_integ.R')
library(UCell)

#Calculating UCell scores
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

get_exp_plot <- function(obj1,obj2,inpath,label='out',
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

obj1 <- readRDS('/data/work/revision/14.revision/03.expression/macaca_pfc_Celltype_finnal.rds')
obj2 <- readRDS('/data/work/revision/14.revision/05.correlation/02.human/human_final_celltype.rds')
obj1$ingroup=obj1$Celltype
obj2$ingroup=obj2$celltype
get_exp_plot(obj1,obj2,label='ma2nhs',inpath='/data/work/revision/14.revision/05.correlation/02.human/')

obj1 <- readRDS('/data/work/revision/14.revision/03.expression/macaca_pfc_Celltype_finnal.rds')
obj2 <- readRDS('/data/work/revision/14.revision/04.macaca_cell/macaca_pfc_SubClass_final.rds')
Idents(obj2) <- 'SubClass'
obj2=subset(obj2,downsample=5000)
obj1$ingroup=obj1$Celltype
obj2$ingroup=obj2$SubClass
get_exp_plot(obj1,obj2,label='ma2nma',inpath='/data/work/revision/14.revision/05.correlation/01.macaca/')

obj1 <- readRDS('/data/work/revision/14.revision/05.correlation/03.hs2hs/human_paper_sub_final.rds')
obj2 <- readRDS('/data/work/revision/14.revision/05.correlation/02.human/human_final_celltype.rds')
obj1$ingroup=obj1$sub
obj2$ingroup=obj2$celltype
get_exp_plot(obj1,obj2,label='hs2nhs',inpath='/data/work/revision/14.revision/05.correlation/03.hs2hs/')

#Comparing coefficient using boxplot
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
get_cor_value <- function(inpath,type='TFs'){
tmp=readRDS(paste0(inpath,'/final_result.rds'))
ob1=tmp[[1]]
ob2=tmp[[2]]
c1=which(ob2 < 0.05)
df=data.frame('value'=ob1[c1])
df$type=type  
return(df)
}

get_mutilcor_value1 <- function(inpath) {
path1=paste0(inpath,'03.tf/')
path2=paste0(inpath,'02.ntf/')
path3=paste0(inpath,'01.ntf/')
df1=get_cor_value(path1,type='TFs')
df2=get_cor_value(path2,type='HVGs')
df3=get_cor_value(path3,type='sHVGs')
df4=Reduce(rbind,list(df1,df2,df3))
return(df4)          
}

inpath='/data/work/revision/14.revision/05.correlation/02.human/'
df=get_mutilcor_value(inpath)
plot_boxplot(df,x='type',y='value', incolor=cor,label='ma2hs')