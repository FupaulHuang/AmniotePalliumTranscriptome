source('functions/function_integ.R')
source("functions/function_merge.R")
source("functions/function_deg.R")
library(gplots)
inpath='/data/work/revision/14.revision/02.degs/02.new_scdeg/'
inlist=c('human','macaca','mouse','zebrafinch','turtle','lizard')
incolor=c("#FED439FF","#709AE1FF","#8A9197FF","#D2AF81FF","#FD7446FF","#D5E4A2FF")
inotg='/data/work/h5/03.graduate_20230712/05.F2/00.ortholog_genes/'

plot_upset_multi(celltype='EX',inlist=inlist)
plot_upset_multi(celltype='IN',inlist=inlist)
plot_upset_multi(celltype='OLI-like lineage',inlist=inlist)
plot_upset_multi(celltype='AST-like_lineage',inlist=inlist)
plot_upset_multi(celltype='Vascular_cell_lineage',inlist=inlist)
plot_upset_mic(celltype='MIC',inlist=c('human','macaca','zebrafinch','turtle','lizard'))
plot_upset_mic(celltype='Vascular_cell_lineage',inlist=c('human','macaca','zebrafinch','turtle','lizard'))

df1=get_upset_multi(celltype='EX',inlist=inlist)
my_venn <- venn(df1, show.plot = FALSE)
#attr(x = my_venn, "intersections")
attr(x = my_venn, "intersections")$`1human:2macaca:3mouse:4zebrafinch:5turtle:6lizard`
attr(x = my_venn, "intersections")$`1human:2macaca:3mouse`
attr(x = my_venn, "intersections")$`4zebrafinch:5turtle:6lizard`