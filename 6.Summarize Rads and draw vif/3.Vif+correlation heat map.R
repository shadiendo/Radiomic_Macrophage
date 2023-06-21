rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(pheatmap)
library(dplyr)
library(corrplot)
library(ggplot2)
library(tidyverse)
library(aplot)

# ████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# VIF
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
df_vif <- read.table(file="→→→VIF_muticoxMiddleRet_num20.csv",
                     header = TRUE, sep = ",")  # ,row.names=1
df_vif$vif <- round(df_vif$vif,2)  # 保留两位小数


# 整理顺序
df_vif$var_name <- factor(df_vif$var_name,
                          levels = rev(df_vif$var_name))


plot_vif <- ggplot(df_vif,aes(var_name,vif))+
  geom_bar(stat = 'identity',position="dodge",fill ="#67001f")+
  geom_text(aes(label=vif),
            position=position_dodge(0.9),  # 当与单个元素宽度不同时的闪避宽度
            # size=4.8,
            vjust=1.2,hjust=0.5,
            color='white'
  )+
  # coord_flip()+  # 翻转轴
  
  
  labs(x=NULL,y='VIF_value')+  # 去掉标签名称
  theme_bw()+
  theme(
    axis.text.x = element_blank(), # 去除刻度标签
    axis.text.y = element_blank(), # 去除刻度标签
    
    axis.ticks = element_blank(),  # 去除刻度线
    panel.border = element_blank(),  # 去掉外边框
    
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank())

plot_vif



# ████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 相关性热图
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

#读取数据
df_orig <- read.table(file="radiomics1402_pretreated_OS.csv",
                      header = TRUE, sep = ",",row.names=1)
feature_selected <- df_vif$var_name


# 根据筛选，提取出原列表中目标列
df_screened <- df_orig[ ,colnames(df_orig) %in% feature_selected]

# 做相关
df_corr = as.data.frame(cor(df_screened))
df_corr[upper.tri(df_corr)]<-NA  # 取一半


# # 拿一下颜色条
# pdf("legendBar.pdf")
df_corr_bar = cor(df_screened)
corrplot(df_corr_bar,
         method = 'color',
         type ='lower',
         tl.pos='n')  # 标签
# dev.off()


# 清理
rm(df_orig,feature_selected,df_screened)

# 变长数据
dfLong <- df_corr %>%
  rownames_to_column("name1") %>%
  pivot_longer(-1,names_to = "name2",values_to = "Value")
dfLong$Value <- round(dfLong$Value,1)  # 保留两位小数
dfLong <- dfLong %>% drop_na(Value)  # 删除所有空值


# 理顺方向
dfLong$name1 <- factor(dfLong$name1,
                       levels = rev(rownames(df_corr))
                       # levels = rownames(df_corr)
)
dfLong$name2 <- factor(dfLong$name2,
                       levels = rev(rownames(df_corr))
)

# ——————————————————————————————————————————————————————————

# 曲线数据
rf_root= paste0(getwd(),'/','前一步的特征/', sep = "")

CoxML="→1000iCoxRet_Selected_num80.txt"
CoxLasso="→1000iCoxLasso_Selected_num35.txt"
Rsfvh="→1000iRSFVH_Selected_num31.txt"

CoxML <- read.table(file=paste0(rf_root,CoxML), sep = ",")
CoxLasso <- read.table(file=paste0(rf_root,CoxLasso), sep = ",")
Rsfvh <- read.table(file=paste0(rf_root,Rsfvh), sep = ",")

# 创建一个曲线数据始祖
n = nrow(df_corr)

# 取对数，减少数据差距
df_curve = df_vif
df_curve$linewidth = log(df_vif$vif+1,base=4)
# 取对称值，把大的变小的，小的变大的
df_curve$linewidth = min(df_curve$linewidth)+max(df_curve$linewidth)-df_curve$linewidth

# 记录每个坐标点的位置
df_curve <- df_curve %>%
  mutate(x1 = c(n:(1)),
         y1 = c(n:1))

# 如果24个特征，那么 21 11 合适
# CoxML 添加目标连接点
CoxML = df_curve %>%
  mutate(x0 = rep(17,n),
         y0 = rep(10,n)) %>%
  filter(var_name %in% CoxML$V1) %>%
  column_to_rownames('var_name')


# 如果24个特征，那么 18 7 合适
# CoxLasso 添加目标连接点
CoxLasso = df_curve %>%
  mutate(x0 = rep(15,n),
         y0 = rep(5,n)) %>%
  filter(var_name %in% CoxLasso$V1) %>%
  column_to_rownames('var_name')

# 如果24个特征，那么 14 3 合适
# Rsfvh 添加目标连接点
Rsfvh = df_curve %>%
  mutate(x0 = rep(10,n),
         y0 = rep(2,n)) %>%
  filter(var_name %in% Rsfvh$V1) %>%
  column_to_rownames('var_name')


# 绘图
plot_coor <- ggplot(dfLong,aes(x=name1,y=name2))+
  geom_tile(aes(fill=Value))+  # 设置热图填充
  geom_text(aes(label=Value))+  # 设置标注字母
  
  # # 设置颜色
  scale_fill_gradientn(colours = rev(c(
    "#67001f",  # -1
    '#b2182b',  # -0.2
    '#d6604d',  # -0.4
    '#f2a280',  # -0.6
    '#fcd9c4',  # -0.8
    "#f9fbfd",  # 0
    '#cfe4ef',  # -0.2
    '#8fc3dd'  # -0.4
    # '#4191c2',  # -0.5
    # '#1f63a8',  # -0.6
    # "#053061"  # 1
  )))+ 
  
  scale_y_discrete(position="left") +      # 设置标注的位置
  
  # 去掉坐标轴名称
  labs(x=NULL,y=NULL)+
  
  theme(
    axis.ticks = element_blank(),  # 去掉小竖线(刻度线)
    panel.background =element_blank(),  # 去掉背景灰色
    
    axis.text.x = element_blank(), # 去除刻度标签
  )+
  
  
  # ——————————————————————————————————————————————————
  geom_curve(data = CoxML,aes(x = x0,y = y0,xend = x1,yend = y1),
             size = CoxML$linewidth,
             color = '#f5b0c8',
             curvature = 0.1)+ #设置弧度
  geom_point(aes(x = CoxML$x0[1],y = CoxML$y0[1]),
             size = 5,color = "#f5b0c8")+
  
  
  geom_curve(data = CoxLasso,aes(x = x0,y = y0,xend = x1,yend = y1),
             size = CoxLasso$linewidth,
             color = '#ec3877',
             curvature = 0.1)+ #设置弧度
  geom_point(aes(x = CoxLasso$x0[1],y = CoxLasso$y0[1]),
             size = 5,color = "#ec3877")+
  
  geom_curve(data = Rsfvh,aes(x = x0,y = y0,xend = x1,yend = y1),
             size = Rsfvh$linewidth,
             color = '#88054d',
             curvature = -0.1)+ #设置弧度
  geom_point(aes(x = Rsfvh$x0[1],y = Rsfvh$y0[1]),
             size = 5,color = "#88054d")+
  
  #添加连线终点的点
  geom_point(data = df_curve,
             aes(x = x1,y = y1),size = 5,color = "black")

plot_coor






plot_ret <- plot_coor %>%
  # insert_right(plot_vif, width = 0.3)
  insert_top(plot_vif, height = 0.2)


plot_ret



ggsave("→3.VifCor_plot.pdf", plot_ret, width = 12, height = 9.3,dpi = 600)#保存文件

