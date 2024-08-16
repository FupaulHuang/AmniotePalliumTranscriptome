
library(Libra)
library(UpSetR)
source("functions/function_merge.R")
source('functions/function_integ.R')
de_intersect_list <- function(fdeg,insplit='method',incol='cluster',fc=1.5,p.val=0.05) {
fdeg=fdeg[fdeg$FC > fc,]
fdeg=fdeg[fdeg$p_val_adj< p.val,]
de_intersect <- function(fdeg,insplit='method',ingroup=NULL){
fdeg=fdeg[fdeg[,'cluster']==ingroup,]
get_list=function(fdeg,intype,insplit='method') {
    lgene=fdeg[,'gene'][fdeg[,insplit]==intype]
    return(lgene)
}
inlist=lapply(unique(fdeg[,insplit]),get_list,fdeg=fdeg)
linter=Reduce(intersect, inlist)
return(linter)
}

lits=lapply(unique(fdeg[,incol]),de_intersect,insplit=insplit,fdeg=fdeg)
names(lits)= unique(fdeg[,incol])
return(lits)
}

sc_de_intersect_list <- function(fdeg,insplit='cluster',fc=1.5,p.val=0.05) {
fdeg=fdeg[fdeg$FC>fc,]
fdeg=fdeg[fdeg$p_val_adj<p.val,]
get_list=function(fdeg,intype,insplit='method') {
    lgene=fdeg[,'gene'][fdeg[,insplit]==intype]
    return(lgene)
}
inlist=lapply(unique(fdeg[,insplit]),get_list,fdeg=fdeg,insplit=insplit)
names(inlist)=unique(fdeg[,insplit])
return(inlist)
}

#https://github.com/hms-dbmi/UpSetR/issues/169
specific_intersections <- function(data, first.col, last.col, intersections, order_mat,
                                   aggregate, decrease, cut, mbar_color, set_names){
  data <- as.data.frame(data)
  sets <- names(data[c(first.col:last.col)])
  keep <- unique(unlist(intersections))
  remove <- sets[which(!sets %in% keep)]
  remove <- which(names(data) %in% remove)
  if(length(remove) != 0){
    data <- data[-remove]
  }

  data <- plyr::count(data[keep])
  sets <- names(data[1:length(keep)])
  data <- lapply(intersections, function(x){
    temp_sets <- unlist(x)
    x <- data[which(rowSums(data[1:length(keep)]) == length(temp_sets)), ]
    x <- x[which(rowSums(x[temp_sets]) == length(temp_sets)), ]
    if(nrow(x) == 0){
      names <- names(x[1:length(keep)])
      x <- rbind(x, rep(0, ncol(x)))
      colnames(x) <- c(names, "freq")
      x[ ,which(names %in% temp_sets)] <- 1
    }
    x <- x
  })
  
  Freqs <- data.frame()
  
  for(i in seq(length(data))){
    Freqs <- rbind(Freqs, data[[i]])
  }
  
  Freqs <- Freqs[c(set_names, "freq")]
  
  num_sets <- length(keep)
  
  if(tolower(aggregate) == "degree" | is.null(order_mat) == TRUE){
    for(i in 1:nrow(Freqs)){
      Freqs$degree[i] <- rowSums(Freqs[ i ,1:num_sets])
    }
    if(is.null(order_mat) == FALSE) { 
        order_cols <- c()
        for(i in 1:length(order_mat)){
          order_cols[i] <- match(order_mat[i], colnames(Freqs))
        }
        
        for(i in 1:length(order_cols)){
          logic <- decrease[i]
          Freqs <- Freqs[order(Freqs[ , order_cols[i]], decreasing = logic), ]
        }
    }
  } else if(tolower(aggregate) == "sets" & is.null(order_mat) == FALSE) {
    Freqs <- Get_aggregates(Freqs, num_sets, order_mat, cut)
  } else {
    stop('Not implemented yet')
  }
  #delete rows used to order data correctly. Not needed to set up bars.
  delete_row <- (num_sets + 2)
  Freqs <- Freqs[ , -delete_row]
  for( i in 1:nrow(Freqs)){
    Freqs$x[i] <- i
    Freqs$color <- mbar_color
  }
  Freqs <- na.omit(Freqs)
  return(Freqs)
}

assignInNamespace('specific_intersections', specific_intersections, ns= 'UpSetR')

plot_upset <- function(celltype='EX'){
ldf=list()
for (i in 1:length(inlist)) {
df1=read.table(paste0(inpath,inlist[i],'_deg_markers.xls'),sep='\t',header=T)
df1=df1[df1$FC>2,]
df1=df1[df1$p_val_adj<0.05,]
df1=subset(df1,cluster==celltype)
if (inlist[i]!='human'){
otg1=read.table(paste0(inotg,'human2',inlist[i],'_orthology_one2one.xls'),sep='\t',header=T)
c1 <- which(df1$gene %in% otg1[,3])
df1=df1[c1,]
c1 <- which(otg1[,3] %in% df1$gene)
otg1=otg1[c1,]
rownames(otg1) <- otg1[,3]
otg1=otg1[df1$gene,]
df1$gene=otg1[,1]
}

ldf[[paste0(i,inlist[i])]]=as.vector(df1[,'gene'])
}
inquery = list(
  list(
    query = intersects,
    params = list('a1human','a2macaca','a3mouse','a4zebrafinch','a5turtle','a6lizard'),
    color = "darkred",
    active = T,
    query.name = "Amniota share"),
  list(
    query = intersects,
    params = list('a1human','a2macaca','a3mouse'),
    color = "darkgreen",
    active = T,
    query.name = "Mammalia share"),
  list(
    query = intersects,
    params = list('a4zebrafinch','a5turtle','a6lizard'),
    active = T,
    color='darkblue',
    query.name = "Sauropsida share")
  )
inset=c('1human','2macaca','3mouse','4zebrafinch','5turtle','6lizard')
indf=fromList(ldf)
colnames(indf) <- paste0('a',colnames(indf))
p1=upset(indf,order.by= NULL,
         sets = c('a1human','a2macaca','a3mouse','a4zebrafinch','a5turtle','a6lizard'),
         intersections=list(
             list('a1human'),
             list('a2macaca'),
             list('a3mouse'),
             list('a1human','a2macaca'),
             list('a1human','a3mouse'),
             list('a2macaca','a3mouse'),
             list('a1human','a2macaca','a3mouse'),
             list('a4zebrafinch'),
             list('a5turtle'),
             list('a6lizard'),
             list('a4zebrafinch','a5turtle'),
             list('a4zebrafinch','a6lizard'),
             list('a5turtle','a6lizard'),
             list('a4zebrafinch','a5turtle','a6lizard'),
             list('a1human','a2macaca','a3mouse','a4zebrafinch','a5turtle','a6lizard')
            ),
         keep.order=T,
         #order.by = c("degree", "degree"),
         point.size = 4,line.size = 2,
         mainbar.y.label =celltype,sets.x.label = "DEGs Per Specie",
         #sort_intersections=FALSE,
         text.scale = c(1.5, 1.5, 1.5, 1.5, 2, 2),
         #decreasing = c(FALSE, FALSE),
         queries = inquery, query.legend = "top",
         main.bar.color='black',
         sets.bar.color=incolor
         #scale.intersections='log2'
         )

#scale.intersections = "log2"
pdf(paste0('upsetplot_',celltype,'.pdf'),10,8)
print(p1)
dev.off()
}

plot_upset_mic <- function(celltype='EX',inlist=NULL){
ldf=list()
inpath1='/data/work/revision/14.revision/02.degs/'
inpath2='/data/work/revision/14.revision/02.degs/02.new_scdeg/'
  
for (i in 1:length(inlist)) {
inspecies=inlist[i]
deg1=read.table(paste0(inpath1,inspecies,'_deg_markers.xls'),sep='\t',header=T)
deg2=read.table(paste0(inpath2,inspecies,'_deg_markers.xls'),sep='\t',header=T)
#AST-like lineage
deg1$cluster=gsub('EG|AST|EPEN','AST-like_lineage',deg1$cluster)
deg2$cluster=gsub('EG|AST|EPEN','AST-like_lineage',deg2$cluster)
#Vascular cell lineage
deg1$cluster=gsub('ENDO|MUR|PERI|VLMC','Vascular_cell_lineage',deg1$cluster)
deg2$cluster=gsub('ENDO|MUR|PERI|VLMC','Vascular_cell_lineage',deg2$cluster)
#OLI-like lineage
deg1$cluster=gsub('OLI|OPC','OLI-like lineage',deg1$cluster)
deg2$cluster=gsub('OLI|OPC','OLI-like lineage',deg2$cluster)

deg1=subset(deg1,method!='edgeR_QLF')
lgene2=sc_de_intersect_list(deg2,insplit='cluster',fc=2,p.val=0.05)
lgene1=de_intersect_list(deg1,insplit='method',incol='cluster',fc=2)
outgene=intersect(lgene1[[celltype]],lgene2[[celltype]])
if (inlist[i]!='human'){
otg1=read.table(paste0(inotg,'human2',inlist[i],'_orthology_one2one.xls'),sep='\t',header=T)
c1 <- which(outgene %in% otg1[,3])
outgene=outgene[c1]
c1 <- which(otg1[,3] %in% outgene)
otg1=otg1[c1,]
rownames(otg1) <- otg1[,3]
otg1=otg1[outgene,]
outgene=otg1[,1]
}

ldf[[paste0(i,inlist[i])]]=as.vector(outgene)
}
inquery = list(
  list(
    query = intersects,
    params = list('a1human','a2macaca','a4zebrafinch','a5turtle','a6lizard'),
    color = "darkred",
    active = T,
    query.name = "Amniota share"),
  list(
    query = intersects,
    params = list('a1human','a2macaca'),
    color = "darkgreen",
    active = T,
    query.name = "Primate share"),
  list(
    query = intersects,
    params = list('a4zebrafinch','a5turtle','a6lizard'),
    active = T,
    color='darkblue',
    query.name = "Sauropsida share")
  )
inset=c('1human','2macaca','4zebrafinch','5turtle','6lizard')
indf=fromList(ldf)
colnames(indf) <- c('a1human','a2macaca','a4zebrafinch','a5turtle','a6lizard')
p1=upset(indf,order.by= NULL,
         sets = c('a1human','a2macaca','a4zebrafinch','a5turtle','a6lizard'),
         intersections=list(
             list('a1human'),
             list('a2macaca'),
             list('a1human','a2macaca'),
             list('a4zebrafinch'),
             list('a5turtle'),
             list('a6lizard'),
             list('a4zebrafinch','a5turtle'),
             list('a4zebrafinch','a6lizard'),
             list('a5turtle','a6lizard'),
             list('a4zebrafinch','a5turtle','a6lizard'),
             list('a1human','a2macaca','a4zebrafinch','a5turtle','a6lizard')
            ),
         keep.order=T,
         #order.by = c("degree", "degree"),
         point.size = 4,line.size = 2,
         mainbar.y.label =celltype,sets.x.label = "DEGs Per Specie",
         #sort_intersections=FALSE,
         text.scale = c(1.5, 1.5, 1.5, 1.5, 2, 2),
         #decreasing = c(FALSE, FALSE),
         queries = inquery, query.legend = "top",
         main.bar.color='black',
         sets.bar.color=incolor[-3]
         #scale.intersections='log2'
         )

#scale.intersections = "log2"
pdf(paste0('upsetplot_',celltype,'.pdf'),8,6)
print(p1)
dev.off()
}

plot_test_deg <- function(celltype='EX',inlist=NULL){
ldf=list()
inpath1='/data/work/revision/14.revision/02.degs/'
inpath2='/data/work/revision/14.revision/02.degs/02.new_scdeg/'
  
for (i in 1:length(inlist)) {
inspecies=inlist[i]
deg1=read.table(paste0(inpath1,inspecies,'_deg_markers.xls'),sep='\t',header=T)
deg2=read.table(paste0(inpath2,inspecies,'_deg_markers.xls'),sep='\t',header=T)
#AST-like lineage
deg1$cluster=gsub('EG|AST|EPEN','AST-like_lineage',deg1$cluster)
deg2$cluster=gsub('EG|AST|EPEN','AST-like_lineage',deg2$cluster)
#Vascular cell lineage
deg1$cluster=gsub('ENDO|MUR|PERI|VLMC','Vascular_cell_lineage',deg1$cluster)
deg2$cluster=gsub('ENDO|MUR|PERI|VLMC','Vascular_cell_lineage',deg2$cluster)
#OLI-like lineage
deg1$cluster=gsub('OLI|OPC','OLI-like lineage',deg1$cluster)
deg2$cluster=gsub('OLI|OPC','OLI-like lineage',deg2$cluster)
deg1=subset(deg1,method!='edgeR_QLF')
lgene2=sc_de_intersect_list(deg2,insplit='cluster',fc=2,p.val=0.05)
lgene1=de_intersect_list(deg1,insplit='method',incol='cluster',fc=2)
outgene=intersect(lgene1[[celltype]],lgene2[[celltype]])
ldf[[paste0(i,inlist[i])]]=length(outgene)/length(lgene2[[celltype]])
}

ldf=t(data.frame(ldf))
return(ldf)
}

run_de_bulk <- function(sc,label='label',cell_type='cell_type',replicate='replicate',incluster='out',FC=1.5){
sc@meta.data[,'label']= sc@meta.data[,label]
sc@meta.data[,'cell_type']= sc@meta.data[,cell_type]
sc@meta.data[,'replicate']= sc@meta.data[,replicate]

run_de_filter <- function(sc,de_method = "edgeR", de_type = "LRT",n_threads = 2,FC=1.5,p_val_adj=0.05,incluster='out') {
de = run_de(sc,de_method = de_method, de_type = de_type, n_threads=n_threads)
de$FC=exp(de$avg_logFC)
de=de[de$FC>FC,]
de=de[de$p_val_adj < p_val_adj,]
de$method=paste0(de$de_method,'_',de$de_type)
de$cluster=incluster
de=data.frame(de)
print(paste0(de_method,'_',de_type,':  ',dim(de)[1]))
return(de)
}
#de2[grep('SLC17A6|SLC17A7',de2$gene),]
n_threads=4
de1 = run_de_filter(sc,de_method = "edgeR",  n_threads=n_threads, incluster=incluster,de_type = "LRT")
de2 = run_de_filter(sc,de_method = "edgeR",  n_threads=n_threads, incluster=incluster,de_type = "QLF")
de3 = run_de_filter(sc,de_method = "DESeq2", n_threads=n_threads, incluster=incluster,de_type = "LRT")
de4 = run_de_filter(sc,de_method = "DESeq2", n_threads=n_threads, incluster=incluster,de_type = "Wald")
de5 = run_de_filter(sc,de_method = "limma",  n_threads=n_threads, incluster=incluster,de_type = "trend")
de6 = run_de_filter(sc,de_method = "limma",  n_threads=n_threads, incluster=incluster,de_type = "voom")

#lde=list(de1$gene,de2$gene,de3$gene,de4$gene,de5$gene,de6$gene)
#de7=Reduce(intersect,lde)
de8=Reduce(rbind,list(de1,de2,de3,de4,de5,de6))

#write.table(de8,'deg.xls',sep='\t',quote=F)
#write.table(de7,'deg_filter.xls',sep='\t',quote=F)
return(de8)
}

run_de_multi <- function(sc,replicate=NULL,cell_type=NULL,incol=NULL,FC=1.5) {
DefaultAssay(sc)='RNA'
set.seed(2)
if (is.null(cell_type)) {
sc$cell_type='test'
}
if (is.null(replicate)) {
sc$replicate='1'
sc$replicate[sample(1:dim(sc)[2], ceiling(dim(sc)[2]/2))]='2'
}

run_de_label<- function(sc,inlabel,incol=NULL,FC=1.5){
sc$label='unknown'
print(paste0('Calculating ',inlabel))
#sc$label[grep('Ex',sc$Celltype)]='EX'
sc$label[grep(inlabel,sc@meta.data[,incol])]=inlabel
sc$label=factor(sc$label,levels=c(inlabel,'unknown'))
deg=run_de_bulk(sc,label='label',cell_type='cell_type',replicate='replicate', incluster=inlabel,FC=FC)
return(deg)
}
ldeg=lapply(unique(sc@meta.data[,incol]),run_de_label,sc=sc,incol=incol,FC=FC)
fdeg=Reduce(rbind,ldeg)
c1=which(is.na(fdeg))
fdeg=fdeg[-c1,]
return(fdeg)
#for (i in unique(sc@meta.data[,incluster]) {}
}

get_upset_multi <- function(celltype='EX',inlist=NULL,inpath1='/data/work/revision/14.revision/02.degs/',
inpath2='/data/work/revision/14.revision/02.degs/02.new_scdeg/'){
ldf=list()
for (i in 1:length(inlist)) {
inspecies=inlist[i]
deg1=read.table(paste0(inpath1,inspecies,'_deg_markers.xls'),sep='\t',header=T)
deg2=read.table(paste0(inpath2,inspecies,'_deg_markers.xls'),sep='\t',header=T)
#AST-like lineage
deg1$cluster=gsub('EG|AST|EPEN','AST-like_lineage',deg1$cluster)
deg2$cluster=gsub('EG|AST|EPEN','AST-like_lineage',deg2$cluster)
#Vascular cell lineage
deg1$cluster=gsub('ENDO|MUR|PERI|VLMC','Vascular_cell_lineage',deg1$cluster)
deg2$cluster=gsub('ENDO|MUR|PERI|VLMC','Vascular_cell_lineage',deg2$cluster)
#OLI-like lineage
deg1$cluster=gsub('OLI|OPC','OLI-like lineage',deg1$cluster)
deg2$cluster=gsub('OLI|OPC','OLI-like lineage',deg2$cluster)

deg1=subset(deg1,method!='edgeR_QLF')
lgene2=sc_de_intersect_list(deg2,insplit='cluster',fc=2,p.val=0.05)
lgene1=de_intersect_list(deg1,insplit='method',incol='cluster',fc=2)
outgene=intersect(lgene1[[celltype]],lgene2[[celltype]])
if (inlist[i]!='human'){
otg1=read.table(paste0(inotg,'human2',inlist[i],'_orthology_one2one.xls'),sep='\t',header=T)
c1 <- which(outgene %in% otg1[,3])
outgene=outgene[c1]
c1 <- which(otg1[,3] %in% outgene)
otg1=otg1[c1,]
rownames(otg1) <- otg1[,3]
otg1=otg1[outgene,]
outgene=otg1[,1]
}

ldf[[paste0(i,inlist[i])]]=as.vector(outgene)
}
return(ldf)
}

plot_upset_multi <- function(celltype='EX',inlist=NULL){
ldf=list()
inpath1='/data/work/revision/14.revision/02.degs/'
inpath2='/data/work/revision/14.revision/02.degs/02.new_scdeg/'
  
for (i in 1:length(inlist)) {
inspecies=inlist[i]
deg1=read.table(paste0(inpath1,inspecies,'_deg_markers.xls'),sep='\t',header=T)
deg2=read.table(paste0(inpath2,inspecies,'_deg_markers.xls'),sep='\t',header=T)
#AST-like lineage
deg1$cluster=gsub('EG|AST|EPEN','AST-like_lineage',deg1$cluster)
deg2$cluster=gsub('EG|AST|EPEN','AST-like_lineage',deg2$cluster)
#Vascular cell lineage
deg1$cluster=gsub('ENDO|MUR|PERI|VLMC','Vascular_cell_lineage',deg1$cluster)
deg2$cluster=gsub('ENDO|MUR|PERI|VLMC','Vascular_cell_lineage',deg2$cluster)
#OLI-like lineage
deg1$cluster=gsub('OLI|OPC','OLI-like lineage',deg1$cluster)
deg2$cluster=gsub('OLI|OPC','OLI-like lineage',deg2$cluster)

deg1=subset(deg1,method!='edgeR_QLF')
lgene2=sc_de_intersect_list(deg2,insplit='cluster',fc=2,p.val=0.05)
lgene1=de_intersect_list(deg1,insplit='method',incol='cluster',fc=2)
outgene=intersect(lgene1[[celltype]],lgene2[[celltype]])
if (inlist[i]!='human'){
otg1=read.table(paste0(inotg,'human2',inlist[i],'_orthology_one2one.xls'),sep='\t',header=T)
c1 <- which(outgene %in% otg1[,3])
outgene=outgene[c1]
c1 <- which(otg1[,3] %in% outgene)
otg1=otg1[c1,]
rownames(otg1) <- otg1[,3]
otg1=otg1[outgene,]
outgene=otg1[,1]
}

ldf[[paste0(i,inlist[i])]]=as.vector(outgene)
}
inquery = list(
  list(
    query = intersects,
    params = list('a1human','a2macaca','a3mouse','a4zebrafinch','a5turtle','a6lizard'),
    color = "darkred",
    active = T,
    query.name = "Amniota share"),
  list(
    query = intersects,
    params = list('a1human','a2macaca','a3mouse'),
    color = "darkgreen",
    active = T,
    query.name = "Mammalia share"),
  list(
    query = intersects,
    params = list('a4zebrafinch','a5turtle','a6lizard'),
    active = T,
    color='darkblue',
    query.name = "Sauropsida share")
  )
inset=c('1human','2macaca','3mouse','4zebrafinch','5turtle','6lizard')
indf=fromList(ldf)
colnames(indf) <- paste0('a',colnames(indf))
p1=upset(indf,order.by= NULL,
         sets = c('a1human','a2macaca','a3mouse','a4zebrafinch','a5turtle','a6lizard'),
         intersections=list(
             list('a1human'),
             list('a2macaca'),
             list('a3mouse'),
             list('a1human','a2macaca'),
             list('a1human','a3mouse'),
             list('a2macaca','a3mouse'),
             list('a1human','a2macaca','a3mouse'),
             list('a4zebrafinch'),
             list('a5turtle'),
             list('a6lizard'),
             list('a4zebrafinch','a5turtle'),
             list('a4zebrafinch','a6lizard'),
             list('a5turtle','a6lizard'),
             list('a4zebrafinch','a5turtle','a6lizard'),
             list('a1human','a2macaca','a3mouse','a4zebrafinch','a5turtle','a6lizard')
            ),
         keep.order=T,
         #order.by = c("degree", "degree"),
         point.size = 4,line.size = 2,
         mainbar.y.label =celltype,sets.x.label = "DEGs Per Specie",
         #sort_intersections=FALSE,
         text.scale = c(1.5, 1.5, 1.5, 1.5, 2, 2),
         #decreasing = c(FALSE, FALSE),
         queries = inquery, query.legend = "top",
         main.bar.color='black',
         sets.bar.color=incolor
         #scale.intersections='log2'
         )

#scale.intersections = "log2"
pdf(paste0('upsetplot_',celltype,'.pdf'),8,6)
print(p1)
dev.off()
}