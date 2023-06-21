rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(reshape)
library(ggplot2)
library(plyr)
library(patchwork)

data <- read.csv("GeneralView - for scalarGram - 肿瘤纯度巨噬细胞.csv")
colnames(data)
# 修改标签名
data <- rename(data, c('Macrophages'='value1','TumorPurity'='value2'))    


#排序
df_order <- data[order(data$RiskScore),] 


# 准备小热图所需的数据
df_cell <- data[,c(1,6)]
# 数据变换,变成三列
df2=melt(df_cell,id="ID") 



#改变level
data$ID <- factor(data$ID, levels = df_order$ID )
df2$ID <- factor(df2$ID, levels = df_order$ID )
df2$variable <- factor(df2$variable,levels = c("RiskScore_log2_Scaling"))



p1 = ggplot(df2,aes(x=ID,y=variable,fill=value))+
  geom_raster()+
  scale_fill_gradient2(high="red3")+   # low="#fadcb1", high="#8b0000", mid="#ed918d"
  scale_y_discrete(position="right") +
  
  # 去除上下两边的空袭
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +

  theme(
    legend.position = "none",
    axis.text.x.bottom = element_blank(),  #修改坐标轴文本大小
    axis.ticks = element_blank(),   #不显示坐标轴刻度
    legend.title = element_blank()   #不显示图例title
  )
p1



# 绘制配对哑铃图并添加拟合线，因为画的是双线，确定下两个变量,分别为
p2 <- ggplot(data) +
  geom_segment(aes(x = ID,xend = ID,y = value1, yend = value2),
               show.legend = TRUE,
               color = "#374e55",    # 连接线的颜色
               size = 0.2,  # 连接线的宽度
               alpha=0.1
               ) +
  xlab(NULL)+
  ylab("Percentage of macrophage infiltration or tumor purity")+

  geom_point(aes(x=ID, y=value1),group = 1,color = "#374e55",size = 3) +
  stat_smooth(aes(x = as.numeric(ID), y = value1),
              method=loess,
              linetype = 2,
              color = '#374e55',
              fill = 'grey',  # 置信区间颜色
              level=0.95) +
  
  geom_point(aes(x=as.numeric(ID), y=value2),color = "#189455",size = 3) +
  stat_smooth(aes(x = as.numeric(ID), y = value2),
              method=loess,
              linetype = 2,
              color = '#189455',
              fill = '#78b99f',   # 置信区间颜色
              level=0.95)+
  theme_classic()+
  theme(
    axis.ticks.x = element_blank(),  # 去掉地下的刷子（刻度线）
    axis.line.x = element_blank(),    # 不知道干啥用
    axis.text.x = element_blank(),  # 去除坐标轴X标签
    axis.ticks = element_blank(),  # 不知道干啥用
    )
p2


# 拼接即可
ret <- p2/p1+plot_layout(heights = c(1, 0.05))

ggsave('肿瘤纯度和巨噬细胞双规散点图.pdf', ret, width = 12, height = 7,dpi = 300)#保存文件











