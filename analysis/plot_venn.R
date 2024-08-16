source("functions/function_plot.R")
library(VennDiagram)
library(dplyr)
library(ggplot2)

df1<- read.table('09.results/xixing_in_deg.xls',header=T)
incolor=c("#440154ff", '#21908dff', '#fde725ff')

df2=df1[grep('low',df1$sp_cluster),]
df3=df1[grep('high',df1$sp_cluster),]
incolor=c("#D5E4A2FF","#FD7446FF","#D2AF81FF")
plot_venn(data=df2,artist='sp_cluster',word='gene',incolor=incolor,label='low')
plot_venn(data=df3,artist='sp_cluster',word='gene',incolor=incolor,label='high')
