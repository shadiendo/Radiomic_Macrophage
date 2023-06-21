rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(pheatmap)
library(dplyr)
library(tidyverse)
library(ggtree)  #聚类
library(aplot)   #拼图
library(patchwork)



# ▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉
# 函数部分
# ▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉

showPlot1<- function(path,titleName){
  # path= '→ALL_p0.05_结果② 分类平均相关度表.csv'
  # titleName = 'DEMO'
  
  # #########################################################
  # 数据整理
  # #########################################################

  # 打开文件内容
  content_orig <- read.table(path, header = TRUE, sep = ",",row.names=1)
  cat('【原始维度】',dim(content_orig),'\n')
  
  # # 取特征
  content <- content_orig[,2:ncol(content_orig)] # 取第4列到最后一列
  cat('【截取维度】',dim(content),'\n')
  
  
  names(content)
  # # 修改列名
  # new_rownames = rbind('RF1','RF2','RF3','RF4','RF5','RF6','RF7',
  #                      'RF8','RF9','RF10','RF11','RF12','RF13',
  #                      'RF14','RF15','RF16','RF17','RF18','RF19','RF20')
  # colnames(content) = new_rownames
  
  
  # 转化为ggplot2格式的数据类型，长数据
  content_long = content %>%
    rownames_to_column("GO_name") %>%
    pivot_longer(-1,names_to = "RFs",values_to = "Value")
  
  
  # 理顺y轴方向，逆向则应该对levels和labels： rev(content_long$GO_ID)
  content_long$GO_name <- factor(content_long$GO_name,
                              levels = content_long$GO_name,
                              labels = content_long$GO_name)
  content_long$RFs <- factor(content_long$RFs,
                             levels = content_long$RFs,
                             labels = content_long$RFs)
  
  # 清理
  rm(path,content_orig)
  
  # #########################################################
  # 基本绘图
  # #########################################################
  
  p = ggplot(content_long,aes(x=RFs,y=GO_name,fill=Value))+
    ggtitle(titleName) +
    
    geom_raster()+
    # scale_fill_gradient2(low="#4575b4", high="#d73027", mid="#fafdc7")+
    
    # scale_fill_gradientn(colours = rev(c(
    #   "#67001f",'#b2182b','#d6604d','#f2a280','#fcd9c4',  # -1 ~ -0.2
    #   "#ffffff",  # 0
    #   '#cfe4ef','#8fc3dd','#4191c2','#1f63a8',"#053061"   # 0.1 ~ 2
    #   ))
    # )+
    
    scale_fill_gradientn(colours = rev(c(
      "#cb1813",'#cb0100','#970000','#630000',
      # '#2c0001',  # -1 ~ -0.2
      "#000101",  # 0
      '#003301','#006502','#009b00','#00c800',"#48dc47"   # 0.1 ~ 2
    ))
    )+
    
    # scale_fill_gradientn(colours = rev(c(
    #   "#a30027",'#d62f29','#f56d43','#ffad61','#fede88',  # -1 ~ -0.2
    #   "#ffffbf",  # 0
    #   '#d8ef89','#a5d869','#63be61','#1c964d',"#006836"   # 0.1 ~ 2
    # ))
    # )+
    
    scale_y_discrete(position="right") +
    theme_minimal()+
    theme(axis.text.x=element_text(angle=90,
                                   hjust = 1,
                                   vjust=1))+
    theme(panel.grid.major=element_blank())
  p

  # #########################################################
  # 聚类树
  # #########################################################
  # complete(默认)、mcquitty、ward.D、ward.D2、single、median、centroid、average
  phY <- content %>%
    dist(method = "manhattan")%>%   # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
    hclust(method = "ward.D2")%>%
    ggtree(layout="rectangular",branch.length="none")
  # phY

  # phX <- content %>%
  #   t()%>%
  #   dist(method = "euclidean")%>%
  #   hclust(method = "complete")%>%
  #   ggtree(layout="rectangular",branch.length="none")+
  #   layout_dendrogram()
  # phX

  p_ret <- p %>%
    # insert_top(pX, height = .05) %>%
    # insert_left(pY, width = .05) %>%
    insert_left(phY,width=.1) 
    # insert_top(phX,height=.1)%>%
  p_ret

  return(p_ret)
  
}






##############################################################
# 平均相关度表
corr_ALL_0.05='→ALL_p0.05_结果② 分类平均相关度表.csv'
corr_ALL_0.05 <- showPlot1(corr_ALL_0.05,'corr_ALL_0.05')


# ——————————————————————————————————————————————————————————————————
ggsave('→corr_ALL_0.05-hotmap.pdf', corr_ALL_0.05, width =10, height = 10, dpi = 300)




