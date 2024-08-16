library(dplyr)
library(ggalluvial)
library(ggrepel)
meta <- read.table('03.cocluster/in_metadata.xls',sep='\t',header=T)
meta$Freq=0.1

meta$first=recode(meta$lab,hs='Human',ma='Macaque',mm='Mouse',zf='Zebra finch',tt='Turtle',li='Lizard')
meta$first=factor(meta$first,levels=rev(c('Human','Macaque','Mouse','Zebra finch','Turtle','Lizard')))
meta$second=paste0('C',meta$seurat_clusters)
meta1=meta
incluster=c('C7','C8','C9')
c1=which(meta1$second %in% incluster)
meta=meta1[c1,]
meta$second=factor(meta$second,levels=paste0('C',rev(c(7:9))))

meta$sub=recode(meta$sub,Sst='SST',Lamp5=toupper('Lamp5'),Pvalb=toupper('Pvalb'),Sncg=toupper('Sncg'),Chodl=toupper('Chodl'),Vip=toupper('Vip'))
meta$sub=gsub('Sst Chodl','SST',meta$sub)
meta$sub=gsub('GABAergic','IN1_tt',meta$sub)

meta1 <- read.table('01.rds/all_metadata_20230714.xls',sep='\t',header=T)
lcell=rownames(meta1)[meta1$Celltype=='tsInh1']
c1 <- which(rownames(meta) %in% lcell)
meta$sub[c1]='IN1_tt'
lcell=rownames(meta1)[meta1$Celltype=='tsInh2']
c1 <- which(rownames(meta) %in% lcell)
meta$sub[c1]='IN2_tt'

meta$sub=recode(meta$sub,InhNeur1='IN1_li',InhNeur2='IN2_li',InhNeur3='IN2_li')

c1 <- which(meta$sub %in% c('GABA-1-1','GABA-1-2'))
meta$sub[c1]='IN1_ZF'
c1 <- which(meta$sub %in% c('GABA-2','GABA-3','GABA-4','GABA-5-1','GABA-5-2','GABA-5-3','GABA-6','GABA-7','GABA-8'))
meta$sub[c1]='IN2_ZF'

meta$sub[grep('IN1',meta$sub)]='Sauropsid_IN1'
meta$sub[grep('IN2',meta$sub)]='Sauropsid_IN2'
meta$sub[grep('PVALB',meta$sub)]='Mammal_PVALB'
meta$sub[grep('VIP|SST|SNCG|PAX6|LAMP5',meta$sub)]='Mammal_IN'

meta$third=factor(meta$sub,levels=rev(unique(meta$sub)))

incolor=c("#FED439FF","#709AE1FF","#8A9197FF","#D2AF81FF","#FD7446FF","#D5E4A2FF")
incolor1=c(incolor,cor[1:14])
incolor=incolor1[-c(7:13)]
ccluster1=c('#8dceba',
           '#ecaf80',
           '#bbb7d8',
           '#f494c4',
           '#b3d28e',
           '#f1d581',
           '#d2ba8e',
           '#b2b2b2',
           '#d2e6f1',
           '#8fbbd8')
ccluster=ccluster1[8:10]
p1=ggplot(data=meta,
       aes(axis1=first,axis2=second,axis3=third,
           y=Freq))+
  geom_alluvium(aes(fill=second),
                #alpha = 0.75,
                #size=3,
                #color="white",
                width = 0.1,
                aes.bind = "flows")+
  scale_fill_manual(values = rev(ccluster))+
  geom_stratum(fill=incolor,
               #color="white",
               #size=3,
               width=0.1)+
  geom_label_repel(stat = "stratum", aes(label = after_stat(stratum)))+
  #geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
  scale_x_continuous(breaks = c(1:3),
                     labels = c("Specie","Cluster",'Cell type'),
                     expand = expansion(mult = c(0,0)))+
  coord_flip()+
  theme_void()

pdf('IN_alluvium_plot_subset1.pdf',8,6)
print(p1)
dev.off()
