library(optparse)
op_list <- list(
make_option(c("-i", "--inrds"), type = "character", default = NULL, action = "store", help = "The input of Seurat RDS",metavar="rds"),
make_option(c("-d", "--ident"), type = "character", default = NULL, action = "store", help = "The sample Ident of Seurat object",metavar="idents"),
make_option(c("-s", "--size"),  type = "integer", default = NULL, action = "store", help = "The sample size of Seurat object",metavar="size"),
make_option(c("-l", "--label"), type = "character", default = "out", action = "store", help = "The label of output file",metavar="label"),
make_option(c("-a", "--assay"), type = "character", default = "Spatial", action = "store", help = "The assay of input file",metavar="assay")
)
parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)

assay <- opt$assay

library(Seurat)
obj <- readRDS(opt$inrds)
if (!is.null(opt$ident)) {
Idents(obj) <-  opt$ident
size=opt$size
if (!is.null(size)) {
obj <- subset(x = obj, downsample = opt$size)
}
saveRDS(obj,"subset.rds")
}
if (is.null(opt$label)) {
label1 <- 'out'
}else{
label1 <- opt$label
}
write.csv(t(as.matrix(obj@assays[[assay]]@counts)),file = paste0(label1,'.csv'),quote=F)
write.table(obj@meta.data,'metadata_subset.xls',sep='\t',quote=F)
