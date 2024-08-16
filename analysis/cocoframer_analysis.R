library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 100000 * 1024^5)

source("functions/function_merge.R")
source('functions/function_integ.R')
source('functions/function_scripts/function_cocoframer.R')
source("functions/function_corplot.R")


obj=readRDS('zf_ex.rds')
obj$third <- obj$ingroup

gene.list=readLines("allTFs_hg38.txt")

#Using TFs
c1 <- which(rownames(obj) %in% gene.list)
obj <- subset(obj,features=rownames(obj)[c1])

DefaultAssay(obj) <- 'RNA'
run_cocoframer(obj,gene.use=NULL,assay.use= 'RNA',group.by='third')

#Using non-TFs/HVGs
c1 <- which(rownames(obj) %in% gene.list)
obj <- subset(obj,features=rownames(obj)[-c1])

DefaultAssay(obj) <- 'RNA'
run_cocoframer(obj,gene.use=NULL,assay.use= 'RNA',group.by='third')

