library(Seurat)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(ggsci)

colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))
len <- 100
cor=c('#a2b88e','#324323')

#turtle
lab <- 'tt_vs_ma-TF'
ob1 <- readRDS('08.corplot/04.ma_vs_tt/03.tf/final_result.rds')
# lab <- 'tt_vs_ma-non_TF'
# ob1 <- readRDS('08.corplot/04.ma_vs_tt/02.ntf/final_result.rds')
ob2 <- ob1[[1]]
v1 <- c('tt_tsInh1_ref','tt_tsInh2_ref')
v2 <- c('ma_VIP_que','ma_SST_que','ma_PVALB_que','ma_LAMP5_que','ma_PAX6_que')
ob2 <- ob2[v2,v1]
g1 <- data.frame(ob2[,1])
g1$group <- 'IN1'
g2 <- data.frame(ob2[,2])
g2$group <- 'IN2'
colnames(g1) <- c('cor_value','group')
colnames(g2) <- c('cor_value','group')
out <- rbind(g1,g2)

my_comparisons1 <- list(c("IN1", "IN2"))
pdf(paste0(lab,'-test.pdf'))
pb1 <-ggboxplot(out,x="group",y='cor_value',color="group",bxp.errorbar.width = 0.5,width = 0.6, size=1)+
      theme_classic() +
      theme(panel.grid.major=element_line(colour=NA),
      plot.title = element_text(hjust = 0.5,size=20),
      axis.title.x=element_text(size=25),
      axis.title.y=element_text(size=25),
      axis.text=element_text(size=24,face = "bold"),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
      legend.text=element_text(colour="black",size=24),
      legend.title=element_text(colour="black",size=25)) +
labs(x = "Celltype", y ="Correlation coefficient", title = paste0(lab)) +
scale_colour_manual(values = cor)
pb1 <- pb1 + stat_compare_means(method="t.test",hide.ns = F,comparisons =my_comparisons1,label="p.format")
print(pb1)
dev.off()

#lizard
lab <- 'li_vs_ma-non_TF'
ob1 <- readRDS('08.corplot/05.ma_vs_li/02.ntf/final_result.rds')
# lab <- 'li_vs_ma-TF'
# ob1 <- readRDS('08.corplot/05.ma_vs_li/03.tf/final_result.rds')
ob2 <- ob1[[1]]
v1 <- c('li_InhNeur1_ref','li_InhNeur2_ref','li_InhNeur3_ref')
v2 <- c('ma_VIP_que','ma_SST_que','ma_PVALB_que','ma_LAMP5_que','ma_PAX6_que')

ob2 <- ob2[v2,v1]
g1 <- data.frame(c(ob2[,1]))
g1$group <- 'IN1'
g2 <- data.frame(c(ob2[,2], ob2[,3]))
g2$group <- 'IN2'
colnames(g1) <- c('cor_value','group')
colnames(g2) <- c('cor_value','group')
out <- rbind(g1,g2)

my_comparisons1 <- list(c("IN1", "IN2"))
pdf(paste0(lab,'-test.pdf'),)
pb1 <-ggboxplot(out,x="group",y='cor_value',color="group",bxp.errorbar.width = 0.5,width = 0.6, size=1)+
      theme_classic() +
      theme(panel.grid.major=element_line(colour=NA),
      plot.title = element_text(hjust = 0.5,size=20),
      axis.title.x=element_text(size=25),
      axis.title.y=element_text(size=25),
      axis.text=element_text(size=24,face = "bold"),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
      legend.text=element_text(colour="black",size=24),
      legend.title=element_text(colour="black",size=25)) +
labs(x = "Celltype", y ="Correlation coefficient", title = paste0(lab)) +
scale_colour_manual(values = cor)


pb1 <- pb1 + stat_compare_means(method="t.test",hide.ns = F,comparisons =my_comparisons1,label="p.format")
print(pb1)

dev.off()


#zebra
lab <- 'li_tt_vs_zf-TF'
ob1 <- readRDS('08.corplot/08.zf_vs_li_tt/03.tf/final_result.rds')
# lab <- 'li_tt_vs_zf-non_TF'
# ob1 <- readRDS('08.corplot/08.zf_vs_li_tt/02.ntf/final_result.rds')
ob2 <- ob1[[1]]

v1 <- c(c('li-tt_InhNeur1_ref','li-tt_InhNeur2_ref','li-tt_InhNeur3_ref'),c('li-tt_tsInh1_ref','li-tt_tsInh2_ref'))
v2 <- rownames(ob2)[grep('_GABA',rownames(ob2))]
v2_1 <- v2[c(4,6)]
v2_2 <- v2[-c(4,6)]

ob2 <- ob2[v2,v1]
g1 <- data.frame(c(ob2[v2_1,1],ob2[v2_1,4]))
g1$group <- 'IN1'
g2 <- data.frame(c(ob2[v2_2,2], ob2[v2_2,3],ob2[v2_2,5]))
g2$group <- 'IN2'
colnames(g1) <- c('cor_value','group')
colnames(g2) <- c('cor_value','group')

out <- rbind(g1,g2)
my_comparisons1 <- list(c("IN1", "IN2"))
pdf(paste0(lab,'-test.pdf'),)
pb1 <-ggboxplot(out,x="group",y='cor_value',color="group",bxp.errorbar.width = 0.5,width = 0.6, size=1)+
      theme_classic() +
      theme(panel.grid.major=element_line(colour=NA),
      plot.title = element_text(hjust = 0.5,size=20),
      axis.title.x=element_text(size=25),
      axis.title.y=element_text(size=25),
      axis.text=element_text(size=24,face = "bold"),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
      legend.text=element_text(colour="black",size=24),
      legend.title=element_text(colour="black",size=25)) +
labs(x = "Celltype", y ="Correlation coefficient", title = paste0(lab)) +
scale_colour_manual(values = cor)


pb1 <- pb1 + stat_compare_means(method="t.test",hide.ns = F,comparisons =my_comparisons1,label="p.format")
print(pb1)
dev.off()

#macaca
lab <- 'ma_vs_zf-non_TF'
ob1 <- readRDS('08.corplot/03.ma_vs_zf/02.ntf/final_result.rds')
# lab <- 'ma_vs_zf-TF'
# ob1 <- readRDS('08.corplot/03.ma_vs_zf/03.tf/final_result.rds')
ob2 <- ob1[[1]]

v1 <- colnames(ob2)[grep('RA_Glut',colnames(ob2))]
v2 <- rownames(ob2)[grep('_L',rownames(ob2))][-8]
v2_1 <- v2[c(2,3,4,5)]
v2_2 <- v2[-c(2,3,4,5)]

ob2 <- ob2[v2,v1]
g1 <- data.frame(c(ob2[v2_1,1],ob2[v2_1,2],ob2[v2_1,3]))
g1$group <- 'upper'
g2 <- data.frame(c(ob2[v2_2,1],ob2[v2_2,2],ob2[v2_2,3]))
g2$group <- 'deep'

colnames(g1) <- c('cor_value','group')
colnames(g2) <- c('cor_value','group')
out <- rbind(g1,g2)

my_comparisons1 <- list(c("upper", "deep"))
pdf(paste0(lab,'-test.pdf'),)
pb1 <-ggboxplot(out,x="group",y='cor_value',color="group",bxp.errorbar.width = 0.5,width = 0.6, size=1)+
      theme_classic() +
      theme(panel.grid.major=element_line(colour=NA),
      plot.title = element_text(hjust = 0.5,size=20),
      axis.title.x=element_text(size=25),
      axis.title.y=element_text(size=25),
      axis.text=element_text(size=24,face = "bold"),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
      legend.text=element_text(colour="black",size=24),
      legend.title=element_text(colour="black",size=25)) +
labs(x = "Celltype", y ="Correlation coefficient", title = paste0(lab)) +
scale_colour_manual(values = cor)

pb1 <- pb1 + stat_compare_means(method="t.test",hide.ns = F,comparisons =my_comparisons1,label="p.format")
print(pb1)
dev.off()

#zebra vs macaca
lab <- 'zf_vs_ma-non_TF'
ob1 <- readRDS('08.corplot/03.ma_vs_zf/02.ntf/final_result.rds')
# lab <- 'zf_vs_ma-TF'
# ob1 <- readRDS('08.corplot/03.ma_vs_zf/03.tf/final_result.rds')
ob2 <- ob1[[1]]

v1 <- colnames(ob2)[grep('_GABA',colnames(ob2))]
v2 <- c('ma_VIP_que','ma_SST_que','ma_PVALB_que','ma_LAMP5_que','ma_PAX6_que')
ob2 <- ob2[v2,v1]
g1 <- data.frame(as.vector(ob2[,c(4,6)]))
g1$group <- 'IN1'
g2 <- data.frame(as.vector(ob2[,-c(4,6)]))
g2$group <- 'IN2'

colnames(g1) <- c('cor_value','group')
colnames(g2) <- c('cor_value','group')

out <- rbind(g1,g2)
my_comparisons1 <- list(c("IN1", "IN2"))

pdf(paste0(lab,'-test.pdf'),)
pb1 <-ggboxplot(out,x="group",y='cor_value',color="group",bxp.errorbar.width = 0.5,width = 0.6, size=1)+
      theme_classic() +
      theme(panel.grid.major=element_line(colour=NA),
      plot.title = element_text(hjust = 0.5,size=20),
      axis.title.x=element_text(size=25),
      axis.title.y=element_text(size=25),
      axis.text=element_text(size=24,face = "bold"),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
      legend.text=element_text(colour="black",size=24),
      legend.title=element_text(colour="black",size=25)) +
labs(x = "Celltype", y ="Correlation coefficient", title = paste0(lab)) +
scale_colour_manual(values = cor)
pb1 <- pb1 + stat_compare_means(method="t.test",hide.ns = F,comparisons =my_comparisons1,label="p.format")
print(pb1)
dev.off()
