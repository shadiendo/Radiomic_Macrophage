rm(list=ls())
# install.packages("reshape2")
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(pheatmap)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggtree)
library(aplot)



titleName='count_ALL'
fileName = '→多维_相关系数表.csv'
data.final <- read.csv(fileName, header = TRUE,check.names=F,row.names = 1)

# 理顺y轴方向，逆向则应该对levels和labels： rev(content_long$GO_ID)
data.final$RFs <- factor(data.final$RFs,
                  levels = data.final$RFs,
                  labels = data.final$RFs)
data.final$ImmuneCell <- factor(data.final$ImmuneCell,
                         levels = data.final$ImmuneCell,
                         labels = data.final$ImmuneCell)


p1 <- ggplot(data.final,aes(x=ImmuneCell,y=RFs))+
  ggtitle(titleName) +
  geom_point(aes(size=`pvalue`,
                 color=`R`))+
  
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,
                                 hjust = 1,
                                 vjust=1))+
  scale_color_gradient(low="#f37704",mid="#f6f5f5",high="#0000cd")+
  # scale_fill_gradient2(low="#1c5ea1", high="#dd725a", mid="#fdfefe")+
  
  
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))
p1


df<-data.final[,c(1,2,3)]   # 聚类用到的是 相关系数R 那一列
# ————————————————————————————————————————————————————
# 以y轴为变量，做层次聚类，并使用ggtree展示层次聚类结果
# ————————————————————————————————————————————————————
# 长格式数据转换为宽格式(热图形式)
df1<-reshape2::dcast(df,
                     RFs ~ ImmuneCell,
                     value.var = "R")
rownames(df1)<-df1$RFs



df1.1<-df1[,2:length(df1)]



# 层次聚类
# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
df1.1.clust<-hclust(dist(df1.1),method='mcquitty')

p2<-ggtree(df1.1.clust)
p2+
  geom_tiplab()+
  xlim(NA,5)


# 使用aplot包拼图
p3 <- p1%>%
  insert_left(p2,width = 0.2)
p3

df1.1.clust$order
df1.1.clust$labels



save_name = paste0 (titleName,'_气泡图.pdf', sep = "")
save_name
# ggsave(save_name, p3, width = 10, height = 7,dpi = 300)#保存文件
