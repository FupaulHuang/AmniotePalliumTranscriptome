source("functions/function_merge.R")
source('functions/function_integ.R')
source('functions/function_plot.R')

plot_idt <- function(df3,label,size=4.5,inlevel,text=F){
df3$specie=factor(df3$specie,levels=inlevel)
df3$variable=factor(df3$variable,levels=inlevel)

theme_heatmap <-
  theme_bw(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(5, "pt"),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(face='bold'),
        axis.text.x = element_text(angle=90,vjust=1, hjust=1))
if (text) {
p1=ggplot(df3,aes(x=specie, y=variable))+
  geom_point(aes(color = value, size = abs(value))) +
  geom_text(aes(label = value), color = 'white', size = 1) +
  scale_size_area(max_size = 5) +
  ggsci::scale_color_material('teal') +
  theme_heatmap +
  theme(legend.key.size = unit(6, "pt"),
        strip.background = element_rect(fill = '#99CCCC'))
    }else{
p1=ggplot(df3,aes(x=specie, y=variable))+
  geom_point(aes(color = value, size = abs(value))) +
  scale_size_area(max_size = 5) +
  ggsci::scale_color_material('teal') +
  theme_heatmap +
  theme(legend.key.size = unit(6, "pt"),
        strip.background = element_rect(fill = '#99CCCC'))    
}
pdf(paste0(label,'_identity.pdf'),size,size)
print(p1)
dev.off()
}

inlevel=c('Human','Cynomolgus monkey','White-tufted-ear marmoset','Mouse',
          'Duckbill platypus',
          'Zebra finch','Rock dove','Chick',
          'Saltwater crocodile','Central bearded dragon','Western painted turtle',
          'African clawed frog','Western clawed frog'
          )
plot_idt(df3,label='STAT6',inlevel=inlevel)

inlevel=c('Human','Cynomolgus monkey','White-tufted-ear marmoset','Mouse',
          'Duckbill platypus',
          'Zebra finch','Rock dove','Chick',
          'Saltwater crocodile','Central bearded dragon','Western painted turtle',
          'African clawed frog','Western clawed frog'
          )
plot_idt(df3,label='NR3C1',inlevel=inlevel)

inlevel=c('Human','Cynomolgus monkey','White-tufted-ear marmoset','Mouse',
          'Duckbill platypus',
          'Zebra finch','Rock dove','Chick',
          'Saltwater crocodile','Central bearded dragon','Western painted turtle',
          'African clawed frog','African clawed frog.1','Western clawed frog'
          )
plot_idt(df3,label='THRB',inlevel=inlevel)
