
library(optparse)
op_list <- list(
make_option(c("-l", "--input_loom"), type = "character", default = NULL, action = "store", help = "The input of aucell loom file",metavar="rds"),
make_option(c("-m", "--input_meta"), type = "character", default = NULL, action = "store", help = "The metadata of Seurat object",metavar="idents"),
make_option(c("-a", "--assay"), type = "character", default = 'Spatial', action = "store", help = "The assay of Seurat object",metavar="assay"),
make_option(c("-c", "--celltype"), type = "character", default = NULL, action = "store", help = "The colname of metadata to calculate RSS",metavar="label")
)
parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)

library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)

celltype <- opt$celltype
assay <- opt$assay
loom <- open_loom(opt$input_loom)

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
close_loom(loom)

meta <- read.table(opt$input_meta,sep='\t',header=T,stringsAsFactor=F)

cellinfo <- meta[,c(opt$celltype,paste0("nFeature_",assay),paste0("nCount_",assay))]
colnames(cellinfo)=c('celltype', 'nGene' ,'nUMI')
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"

sub_regulonAUC <- regulonAUC
rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])
rss=na.omit(rss)
try({
rssPlot <- plotRSS(rss)
save(regulonAUC,rssPlot,regulons,file='regulon_RSS.Rdata')
})

saveRDS(rss,paste0(celltype,"_rss.rds"))

source('00.scripts/function_pyscenic_visualize.R')
plot_pyscenic(inloom='aucell.loom',incolor=incolor,inrss=paste0(celltype,"_rss.rds"),inrds='subset.rds',infun='median', ct.col=celltype,inregulons=NULL,ingrn='grn.tsv',ntop1=5,ntop2=50)
