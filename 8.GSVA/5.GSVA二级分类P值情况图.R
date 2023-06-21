rm(list=ls())
# install.packages('patchwork')
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(pheatmap)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)



showPlot2<- function(fileName,titleName){
  
  # ###############################################
  # 读取数据
  # ###############################################
  
  root= paste0 (getwd(),'/', sep = "", collapse = NULL)
  
  # fileName='TCGA/GSVA_TCGA_二维_相关系数表_p0.05_最终瘦子表.csv'
  # titleName='count_TCGA_0.05'

  # 读取文件
  path= paste0 (root,fileName, sep = "", collapse = NULL)
  cat('【处理文件】:\n',path,'\n')
  
  data.final <- read.csv(path, header = TRUE,check.names=F)
  cat(colnames(data.final))

  # 理顺y轴方向
  data.final$ID <- factor(data.final$ID,
                                 levels = data.final$ID,
                                 labels = data.final$ID)
  
  rm(root,fileName,path)
  # ###############################################
  # 绘图
  # ###############################################

  plot_pop <- ggplot(data.final,aes(x=class,y=ID))+
    ggtitle(titleName) +
    geom_point(aes(size=`Count_log`,
                   color=`coor`))+
    
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x=element_text(angle=0,
                                   hjust = 1,
                                   vjust=0.5))+
    scale_color_gradient(low="lightgrey",high="blue")+
    labs(x=NULL,y=NULL)+
    guides(size=guide_legend(order=3))
  
  plot_pop
  
  # save_name = paste0 (titleName,'_气泡图.pdf', sep = "")
  # save_name
  # ggsave(save_name, plot_pop, width = 7, height = 13,dpi = 300)#保存文件
  

  cat('绘制完毕')
  
  return(plot_pop)
  
}





##############################################################
# 瘦表
count_TCGA_0.05 = '→GSVA后数据/TCGA/GSVA_TCGA_二维_相关系数表_p0.05_最终瘦子表.csv'
count_TCGA_lgg_0.05 = '→GSVA后数据/TCGA_LGG/GSVA_TCGA_LGG_二维_相关系数表_p0.05_最终瘦子表.csv'
count_TCGA_gbm_0.05 = '→GSVA后数据/TCGA_GBM/GSVA_TCGA_GBM_二维_相关系数表_p0.05_最终瘦子表.csv'

count_RBRT_0.05 = '→GSVA后数据/RBRT/GSVA_RBRT_二维_相关系数表_p0.05_最终瘦子表.csv'
count_RBRT_lgg_0.05 = '→GSVA后数据/RBRT_LGG/GSVA_RBRT_LGG_二维_相关系数表_p0.05_最终瘦子表.csv'
count_RBRT_gbm_0.05 = '→GSVA后数据/RBRT_GBM/GSVA_RBRT_GBM_二维_相关系数表_p0.05_最终瘦子表.csv'

count_JSPH_0.05='→GSVA后数据/JSPH/GSVA_JSPH_二维_相关系数表_p0.05_最终瘦子表.csv'

count_ALL_0.05='→GSVA后数据/ALL/GSVA_ALL_二维_相关系数表_p0.05_最终瘦子表.csv'
count_ALL_lgg_0.05='→GSVA后数据/ALL_LGG/GSVA_ALL_LGG_二维_相关系数表_p0.05_最终瘦子表.csv'
count_ALL_gbm_0.05='→GSVA后数据/ALL_GBM/GSVA_ALL_GBM_二维_相关系数表_p0.05_最终瘦子表.csv'

##############################################################
# 瘦表
count_TCGA_0.05 = showPlot2(count_TCGA_0.05,'count_TCGA_0.05')
count_TCGA_lgg_0.05 = showPlot2(count_TCGA_lgg_0.05,'count_TCGA_lgg_0.05')
count_TCGA_gbm_0.05 = showPlot2(count_TCGA_gbm_0.05,'count_TCGA_gbm_0.05')

count_RBRT_0.05 = showPlot2(count_RBRT_0.05,'count_RBRT_0.05')
count_RBRT_lgg_0.05 = showPlot2(count_RBRT_lgg_0.05,'count_RBRT_lgg_0.05')
count_RBRT_gbm_0.05 = showPlot2(count_RBRT_gbm_0.05,'count_RBRT_gbm_0.05')

count_JSPH_0.05=showPlot2(count_JSPH_0.05,'count_JSPH_0.05')

count_ALL_0.05=showPlot2(count_ALL_0.05,'count_ALL_0.05')
count_ALL_lgg_0.05=showPlot2(count_ALL_lgg_0.05,'count_ALL_lgg_0.05')
count_ALL_gbm_0.05=showPlot2(count_ALL_gbm_0.05,'count_ALL_gbm_0.05')


demo <- (count_TCGA_0.05 + count_TCGA_lgg_0.05 + count_TCGA_gbm_0.05) /
  (count_RBRT_0.05+count_RBRT_lgg_0.05+count_RBRT_gbm_0.05) /
  (count_ALL_0.05+count_ALL_lgg_0.05+count_ALL_gbm_0.05)

ggsave('气泡图.pdf', demo, width =20, height = 32, dpi = 600)






