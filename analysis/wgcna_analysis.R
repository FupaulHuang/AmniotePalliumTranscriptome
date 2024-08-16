library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 100000 * 1024^5)
args <- commandArgs(T)
input <- args[1]

source('functions/function_wgcna.R')

obj <- readRDS('input.rds')
seurat_wgcna(obj,group.by='sub',size=100,slot='data',assay='RNA',fun=mean,
                         type = "unsigned",corType = "pearson", pct.mad=0.75, label='test',trait=NULL,nThreads=10)

ul=get_module(inpath='01.upper',label='UL',ingene=neo2dc)
#THRB STAT6
dl=get_module(inpath='02.deeper',label='DL',ingene=neo2dc)
#NR3C1 STAT6
dc=get_module(inpath='05.dc',label='DC',ingene=neo2dc)
#NR3C1 SREBF2 MAFK ARID3A SREBF1 BACH2 GABPA NFE2L1 USF1 STAT6 THRB EP300 TBP ATF4 ELF2 MAFF
c1=which(dc$Node1 %in% c('NR3C1','BACH2','NFE2L1','STAT6','THRB','ATF4','ELF2'))
c2=which(dc$Module1 %in% unique(dc$Module1[c1]))
dc=dc[c2,]

gene1=c('BDNF','GFRA1')
gene2=c('COL6A3','ADAMTS18','NFATC1')
dvr=get_module(inpath='06.dvr',label='DVR',ingene=neo2dvr)
#BCLAF1
hvc=get_module(inpath='03.hvc',label='HVC',ingene=gene1)
#GFRA1
ra=get_module(inpath='04.ra',label='RA',ingene=neo2dvr)
#ELK4

mdc=Reduce(rbind,list(ul,dl,dc))
mdvr=Reduce(rbind,list(dvr,hvc,ra))

#plot module
library(Seurat)
all <- readRDS('01.EX/6libs_integ.rds')
meta=read.table('11.ex_plot/ex_metadata_20231201.xls',sep='\t',header=T)
meta$type=meta$third
meta$type=gsub('LC','others',meta$type)
c1 <- which(meta$type %in% c('L2_IT','L2/3_IT','L3/4_IT','L4_IT'))
meta$type[c1]='UL'
meta$type[grep('^L',meta$type)]='DL'
c1 <- which(meta$type %in% c('AMY','CLA-like','HIP-like'))
meta$type[c1]='others'
all@meta.data=meta
all=subset(all,subset=type!='others')
all@meta.data$type=factor(all@meta.data$type,levels=c('UL','DL','DC','DVR','HVC','RA'))

get_score(all,inlist=mdc,label='DC')
get_score(all,inlist=mdvr,label='DVR')

df1=read.table('DC_addscore.xls',sep='\t',header=T)
plot_module(df1=df1,group_by='first',label='DC')
plot_module(df1=df1,group_by='second',label='DC')
plot_module(df1=df1,group_by='third',label='DC')
plot_module(df1=df1,group_by='type',label='DC')

df1=read.table('DVR_addscore.xls',sep='\t',header=T)
plot_module(df1=df1,group_by='first',label='DVR')
plot_module(df1=df1,group_by='second',label='DVR')
plot_module(df1=df1,group_by='third',label='DVR')
plot_module(df1=df1,group_by='type',label='DVR')
