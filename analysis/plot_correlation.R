library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 100000 * 1024^5)


source("functions/function_corplot.R")
obj1 <- readRDS('/data/work/revision/14.revision/03.expression/macaca_pfc_Celltype_finnal.rds')
obj2 <- readRDS('/data/work/revision/14.revision/04.macaca_cell/macaca_pfc_SubClass_final.rds')
Idents(obj2) <- 'SubClass'
obj2=subset(obj2,downsample=5000)
obj1$ingroup=obj1$Celltype
obj2$ingroup=obj2$SubClass

col1 <- 'ingroup'
col2 <- 'ingroup'
lab1 <- 'ma'
lab2 <- 'tma'

#Calculating correlations using TFs
hvg=NULL
invert=NULL
gene.list <- readLines('/data/work/graduation//09.summary/00.homology_genes/allTFs_hg38.txt')
run_corplot(ob1=obj1,ob3=obj2,ig1=NULL,ig2=NULL,col1=col1,col2=col2,lab1=lab1,lab2=lab2,gene.list=gene.list,assay1='RNA',assay2='RNA',deg=NULL,hvg=hvg,invert=invert)

#Calculating correlations using HVGs
hvg=3000
invert=T
gene.list <- readLines('/data/work/graduation//09.summary/00.homology_genes/allTFs_hg38.txt')
run_corplot(ob1=obj1,ob3=obj2,ig1=NULL,ig2=NULL,col1=col1,col2=col2,lab1=lab1,lab2=lab2,gene.list=gene.list,assay1='RNA',assay2='RNA',deg=NULL,hvg=hvg,invert=invert)

#Calculating correlations using sHVGs
hvg=5000
invert=T
gene.list <- readLines('/data/work/graduation//09.summary/00.homology_genes/allTFs_hg38.txt')
run_corplot(ob1=obj1,ob3=obj2,ig1=NULL,ig2=NULL,col1=col1,col2=col2,lab1=lab1,lab2=lab2,gene.list=gene.list,assay1='RNA',assay2='RNA',deg=NULL,hvg=hvg,invert=invert,shvg=T)