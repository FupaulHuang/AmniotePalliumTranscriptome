
# **********************************************************
# * Author        : HuangFubaoqian
# * Email         : huangbaoqian@genomics.cn
# * Create time   : 2022-10-31 14:36
# * Filename      : plot_featureplot.R
# * Description   : 
# **********************************************************
.libPaths(c("/jdfssz1/ST_SUPERCELLS/P22Z10200N0639/huangbaoqian/software/R4.0/lib"))
source("/jdfssz3/ST_STOMICS/P22Z10200N0664/2.CrossSpecies_Spatial/huangbaoqian/02.data/03.20220504/03.cluster/function_merge.R")
source('/jdfssz3/ST_STOMICS/P22Z10200N0664/2.CrossSpecies_Spatial/huangbaoqian/02.data/03.20220504/03.cluster/04.integ/function_integ.R')


library(optparse)
op_list <- list(
make_option(c("-i", "--inrds"), type = "character", default = NULL, action = "store", help = "The input of Seurat RDS",metavar="rds"),
make_option(c("-d", "--genelist"), type = "character", default = NULL, action = "store", help = "The genelist to plot",metavar="idents"),
make_option(c("-s", "--size"),  type = "integer", default = 1.5, action = "store", help = "The point size ",metavar="size"),
make_option(c("-l", "--label"), type = "character", default = "out", action = "store", help = "The label of output file",metavar="label"),
make_option(c("-a", "--assay"), type = "character", default = "Spatial", action = "store", help = "The assay of input file",metavar="assay"),
make_option(c("-b", "--subset"), type = "character", default = NULL, action = "store", help = "subset rds",metavar="subset")
)
parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)
dir.create(paste0('out_',opt$label))

all <- readRDS(opt$inrds)

n1 <- dim(all)[2]
if (!is.null(opt$subset)) {
library(stringr)
sub <- opt$subset
sub1 <- unlist(str_split(sub,'-'))
sub <- sub1[-1]
sub1 <- sub1[1]

c1 <- which(all@meta.data[,sub1] %in% sub)
all <- subset(all,cells=Cells(all)[c1])
}

n2 <- dim(all)[2]

DefaultAssay(all) <- opt$assay

blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)

fp_plot <- function(obj,gene.list=NULL,output="test",raster=NULL,color=blu,pt.size=1) {
all <- obj
gene <- unique(gene.list)
c <- which(gene %in% rownames(obj))
gene <- gene[c]
for (i in 1:length(gene)){
pf <- paste0(i,"_",gene[i],"_Featureplot.pdf")
pf <- gsub("/","_",pf)

cal_wvsh <- function(obj){
wid <- max(obj@reductions$umap@cell.embeddings[,1])-min(obj@reductions$umap@cell.embeddings[,1])
hei <- max(obj@reductions$umap@cell.embeddings[,2])-min(obj@reductions$umap@cell.embeddings[,2])
rat <- wid/hei
return(rat)
}
rat <- cal_wvsh(all)
wei=rat*10
hei=10

pdf(paste0(output,'_',pf),wei,hei)
if (is.null(raster)) {
print(FeaturePlot(obj,features = gene[i],pt.size=pt.size)+scale_color_gradientn(colours = color)+NoAxes())
}else{
print(FeaturePlot(obj,features = gene[i], raster=F,pt.size=pt.size)+scale_color_gradientn(colours = color)+NoAxes())
}
if (is.null(raster)) {
print(FeaturePlot(obj,features = gene[i],pt.size=pt.size)+scale_color_gradientn(colours = color)+NoAxes()+NoLegend())
}else{
print(FeaturePlot(obj,features = gene[i], raster=F,pt.size=pt.size)+scale_color_gradientn(colours = color)+NoAxes()+NoLegend())
}

dev.off()
}

}

lgene <- readLines(opt$genelist)

psize <- as.numeric(opt$size)

#if (n1!=n2) {
if (FALSE) {
psize=(12+log(1/n2))/2+1
#psize=log(n1)+log(1/n2)+1
print(psize)
}
setwd(paste0('out_',opt$label))
fp_plot(all,gene.list=lgene,output=opt$label,color=blu[-c(1:10)],pt.size=psize)
