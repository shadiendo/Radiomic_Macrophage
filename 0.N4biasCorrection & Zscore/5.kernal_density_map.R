# █████████████████████████████████████████████████████████████████████████████

# 核密度图

# █████████████████████████████████████████████████████████████████████████████
rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1.数据准备
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
df_N4 <- read.csv('→直方图数据-N4.csv', header=T, row.names=1)
df_N4_Zscored <- read.csv('→直方图数据-N4-Zscored.csv', header=T, row.names=1)
df_N4_WSed <- read.csv('→直方图数据-N4-WSed.csv', header=T, row.names=1)

# colnames(df_orig)
# fix(df_orig)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2.绘图函数
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
showPlot<- function(df_demo,titleName_demo,haveLegend,X_limit_min,X_limit_max){
  # df_demo = df_orig
  # titleName_demo='title'
  # haveLegend = 'closeLegend'
  
  
  # 转化为ggplot2格式的数据类型，长数据
  df_demo_long = df_demo %>%
    rownames_to_column("RF_names") %>%
    pivot_longer(-1,names_to = "class",values_to = "Value")
  
  # # # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # # # 绘制核密度图 ggplot 单线
  # # # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # ggplot(df_orig) +
  #   geom_density(aes(segmentation.48,y=..density..),size=1,color="black")+
  #   # theme_minimal()+
  #   theme(panel.grid.major=element_blank())+
  #   xlim(-100, 100)
  
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # 绘制核密度图 ggplot 多线
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  gg <- ggplot(df_demo_long) +
    ggtitle(titleName_demo)+
    xlab("Number Of Repeats") +
    ylab("Kernel Density") +
    geom_density(alpha=.1,
                 size = 0.1,
                 aes(x=Value,
                     # color = class,  # 线条颜色
                     fill = class    # 填充的颜色
                 )
    )+
    xlim(X_limit_min,X_limit_max)+
    # xlim(-30, 50)+
    # coord_cartesian(ylim = c(0,0.05),expand = TRUE)+     # 缩放y轴局部
    
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    theme(plot.title = element_text(color="black", size=18, face="bold",hjust=0.5),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.y = element_text(size = 12,colour="black"),
          axis.text.x = element_text(size = 12,colour="black"))
  
  if (haveLegend=='closeLegend'){
    gg = gg+scale_fill_discrete(guide=FALSE)
  }else if (haveLegend=='openLegend'){
    gg = gg+theme(legend.position = "bottom")
    cat('show legend')
  }
  return(gg)
  
}


p1 = showPlot(df_N4,'N4','closeLegend',0,2000)
p2 = showPlot(df_N4_Zscored,'N4_Zscored','closeLegend',-0.5,4.5)
p3 = showPlot(df_N4_WSed,'N4_Wsed','closeLegend',-75,100)

ret = p1|p2|p3
# ret



ggsave(limitsize = FALSE,'→核密度图结果.pdf',
  ret, width =15, height = 5, dpi = 300)



# 
# pnnn <- showPlot(df_N4_Zscored,'df_N4_WSed','openLegend',-2,6)
# ggsave(limitsize = FALSE,'→中间.pdf',
       # pnnn, width =20, height = 10, dpi = 300)





