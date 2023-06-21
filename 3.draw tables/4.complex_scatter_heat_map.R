rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(reshape)
library(ggplot2)
library(plyr)
library(patchwork)

data <- read.csv("GeneralView - for scalarGram.csv")
data <- data[data$class == 'test', ]  # 只要训练集数据  train   test
colnames(data)

# 以RiskScore排序
df_order <- data[order(data$RiskScore),] 
# 调整level
data$ID <- factor(data$ID, levels = df_order$ID )



# ████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 准备RiskScore条所需的数据
df_cell_radscore <- data[,c("ID","RiskScore_ord")]
df_radscore_order <- df_cell_radscore[order(df_cell_radscore$RiskScore_ord),] 
# 调整level
df_cell_radscore$ID <- factor(df_cell_radscore$ID, levels = df_radscore_order$ID )

# 数据变换,变成三列
df2_radscore=melt(df_cell_radscore,id="ID") 

#改变level
df2_radscore$variable <- factor(df2_radscore$variable,levels = c("RiskScore_ord"))

p0 <- ggplot(df2_radscore,aes(x=ID,y=variable,fill=value))+
  geom_raster()+
  scale_fill_gradient2(low="#fadcb1", mid="#fadcb1",high="red3")+   # low="#fadcb1", high="#8b0000", mid="#ed918d"
  scale_y_discrete(position="right") +
  ylab(NULL)+
  xlab(NULL)+
  
  # 去除上下两边的空袭
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  
  theme(
    legend.position = "none",
    axis.text.x.bottom = element_blank(),  #修改坐标轴文本大小
    axis.ticks = element_blank(),   #不显示坐标轴刻度
    legend.title = element_blank()   #不显示图例title
  )
# p0



# ████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# # 准备小热图所需的数据
# vvvv_needed = c("class_big",
#                 "gender","age",
#                 "Histology","Grade",
#                 "IDH","co_del_1p_19q","TERT_promoter_status","MGMT_status",
#                 "Macrophages","TumorPurity",
#                 "LR","centerline","location_frontal","location_temproal","location_Insular",
#                 "location_parietal","location_occipital","location_other")
# df_cell <- data[,c("ID",vvvv_needed)]
# 
# # 数据变换,变成三列的长数据
# df2=melt(df_cell,id="ID") 
# 
# # 调整level
# df2$ID <- factor(df2$ID, levels = df_order$ID )
# df2$variable <- factor(df2$variable,levels = rev(vvvv_needed))
# 
# #设置颜色
# cols=c(
#   'JSPH'='#80b49b','TCGA_LGG'='#3688be','TCGA_GBM'='#093370','RBRT'='#f0f0f0',
#   # 基线
#   'male'='#000000','female'='#666666',
#   'Q1_age'='#cccccc','Q2_age'='#666666','Q3_age'='#3c3c3c','Q4_age'='#000000',
#   # 病理
#   # 'Astrocytoma'='#efce88','Oligodendroglioma'='#d77d5b','Oligoastrocytoma'='#bc3804','Glioblastoma'='#662800',
#   # 'II'='#feeea9','III'='#9a3d09','IV'='#662800',
#   'Astrocytoma'='#5b5bdd','Oligodendroglioma'='#aacfe4','Oligoastrocytoma'='#3688be','Glioblastoma'='#093370',
#   'II'='#aacfe4','III'='#3688be','IV'='#093370',
#   # 分子病理
#   'Mutant'='#000000','Wild-type'='#666666',
#   'Cod'='#000000','Non-Cod'='#666666',
#   'Methylated'='#000000','unMethylated'='#666666',
#   # 计算结果
#   'Q1_m'="#f3f9fe",'Q2_m'="#aacfe4",'Q3_m'="#3688be",'Q4_m'="#093370",
#   'Q1_tp'="#f3f9fe",'Q2_tp'="#aacfe4",'Q3_tp'="#3688be",'Q4_tp'="#093370",
#   # 位置
#   'R'='#3c3c3c','L'='#666666','LR'='#000000',
#   'YES'='#000000','no'='#666666',
#   
#   # 一般
#   'Not_Available'='#f0f0f0','NOS'='#f0f0f0','--'='#f0f0f0'
#   )
# 
# # 最好关掉geom_tile中的参数
# p1 <- ggplot(df2,aes(x=ID,y=variable,fill=value))+
#   geom_tile(
#     # aes(fill=value),
#     # color='#f0f0f0',
#     # linewidth=F
#   )+
#   xlab(NULL)+
#   ylab(NULL)+
#   
#   scale_x_discrete(expand = c(0,0)) +
#   scale_y_discrete(expand = c(0,0)) +
#   
#   scale_fill_manual(values = cols)+ #指定自定义的颜色
#   
#   theme(
#     legend.position = "bottom",
#     axis.text.x.bottom = element_blank(),  #修改坐标轴文本大小
#     axis.ticks = element_blank(),   #不显示坐标轴刻度
#     legend.title = element_blank(),   #不显示图例title
#     panel.spacing.y = unit(1, "mm")   # 不知道有什么用
#   )
# p1

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
plot_manually <- function(df_demo,Y_title,nrow_demo,ncol_demo) {
  plot_ret <<- ggplot(df_demo,aes(x=ID,y=variable,fill=value))+
    geom_tile(
      # aes(fill=value),
      # color='#f0f0f0',
      # linewidth=F
    )+
    xlab(NULL)+
    ylab(NULL)+
    labs(fill = Y_title)+
    
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    
    scale_fill_manual(values = cols)+ #指定自定义的颜色
  
    guides(fill = guide_legend(nrow = nrow_demo, ncol = ncol_demo))+
    theme(
      axis.text.x.bottom = element_blank(),  #修改坐标轴文本大小
      axis.ticks = element_blank(),   #不显示坐标轴刻度
      panel.spacing.y = unit(0, "cm"),   # 绘图区域中相邻面板之间的垂直间距
      
      legend.position = "right",
      legend.direction = 'vertical',  # horizontal，vertical
      legend.key.width = unit(0.3, "cm"),
      legend.key.height = unit(0.3, "cm"),
      # legend.title = element_blank(),   #不显示图例title
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      # legend.margin = margin(t = 0, r = 1, b = 0, l = 0, unit = "cm")
    )
  return(plot_ret)
}



# 设置基线
vvvv_needed = c("class_big")
df_cell <- data[,c("ID",vvvv_needed)]

# 数据变换,变成三列的长数据
df2=melt(df_cell,id="ID")

# 调整level
df2$ID <- factor(df2$ID, levels = df_order$ID )
df2$variable <- factor(df2$variable,levels = rev(vvvv_needed))

#设置颜色
cols=c(
  'JSPH'='#80b49b','TCGA_LGG'='#3688be','TCGA_GBM'='#093370','RBRT'='#f0f0f0')

p_group <- plot_manually(df2,'Dataset',1,4)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 设置基线
vvvv_needed = c("gender","age")
df_cell <- data[,c("ID",vvvv_needed)]

# 数据变换,变成三列的长数据
df2=melt(df_cell,id="ID")

# 调整level
df2$ID <- factor(df2$ID, levels = df_order$ID )
df2$variable <- factor(df2$variable,levels = rev(vvvv_needed))

#设置颜色
cols=c(
  'JSPH'='#80b49b','TCGA_LGG'='#3688be','TCGA_GBM'='#093370','RBRT'='#f0f0f0',
  # 基线
  'male'='#000000','female'='#666666',
  'Q1_age'='#cccccc','Q2_age'='#666666','Q3_age'='#3c3c3c','Q4_age'='#000000',
  # 一般
  'Not_Available'='#f0f0f0','NOS'='#f0f0f0','--'='#f0f0f0')

p_baseline <- plot_manually(df2,'Baseline',2,10)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 设置病理
vvvv_needed = c("Histology","Grade",
                "IDH","co_del_1p_19q")
df_cell <- data[,c("ID",vvvv_needed)]

# 数据变换,变成三列的长数据
df2=melt(df_cell,id="ID")

# 调整level
df2$ID <- factor(df2$ID, levels = df_order$ID )
df2$variable <- factor(df2$variable,levels = rev(vvvv_needed))

#设置颜色
cols=c(
  # 'Astrocytoma'='#efce88','Oligodendroglioma'='#d77d5b','Oligoastrocytoma'='#bc3804','Glioblastoma'='#662800',
  # 'II'='#feeea9','III'='#9a3d09','IV'='#662800',
  'Astrocytoma'='#5b5bdd','Oligodendroglioma'='#aacfe4','Oligoastrocytoma'='#3688be','Glioblastoma'='#093370',
  'II'='#aacfe4','III'='#3688be','IV'='#093370',
  # 分子病理
  'Wild-type'='#000000','Mutant'='#666666',
  'Cod'='#000000','Non-Cod'='#666666',
  'Methylated'='#000000','unMethylated'='#666666',
  # 一般
  'Not_Available'='#f0f0f0','NOS'='#f0f0f0','--'='#f0f0f0'
  )

p_pathology <- plot_manually(df2,'Pathology',4,10)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 设置位置
vvvv_needed = c("LR","centerline","location_frontal","location_temproal","location_Insular",
                "location_parietal","location_occipital","location_other")
df_cell <- data[,c("ID",vvvv_needed)]

# 数据变换,变成三列的长数据
df2=melt(df_cell,id="ID") 

# 调整level
df2$ID <- factor(df2$ID, levels = df_order$ID )
df2$variable <- factor(df2$variable,levels = rev(vvvv_needed))

#设置颜色
cols=c(
  'R'='#3c3c3c','L'='#666666','LR'='#000000',
  'YES'='#3c3c3c','no'='#666666',
  
  # 一般
  'Not_Available'='#f0f0f0','NOS'='#f0f0f0','--'='#f0f0f0'
)

p_location <- plot_manually(df2,'Region',3,10)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 设置RNAseq
vvvv_needed = c("Macrophages","TumorPurity")
df_cell <- data[,c("ID",vvvv_needed)]

# 数据变换,变成三列的长数据
df2=melt(df_cell,id="ID")

# 调整level
df2$ID <- factor(df2$ID, levels = df_order$ID )
df2$variable <- factor(df2$variable,levels = rev(vvvv_needed))

#设置颜色
cols=c(
  'Q1_m'="#aacfe4",'Q2_m'="#3688be",'Q3_m'="#093370",'Q4_m'="#000000",
  'Q1_tp'="#aacfe4",'Q2_tp'="#3688be",'Q3_tp'="#093370",'Q4_tp'="#000000",
  # 一般
  'Not_Available'='#f0f0f0','NOS'='#f0f0f0','--'='#f0f0f0'
  )
#80b49b
p_rna <- plot_manually(df2,'bioinformatics analysis',2,10)

# ████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 绘制散点图并添加拟合线，因为画的是双线，确定下两个变量,分别为
p_Scatter <- ggplot(data,aes(x=RiskScore_ord, y=time,color=status)) +
  geom_point(size=1.5,stroke = 0)+
  # scale_color_manual(values = c("#99b96d", "#5592c9"))+
  stat_smooth(method=loess) +    # loess 表示线性回归
  
  coord_cartesian(ylim = c(1, 60))+
  scale_x_discrete(expand = c(0.01, 0))+
  scale_y_discrete(limits = c(1,6,12,24, 36, 60),
                   expand = c(0.01, 0))+
  xlab(NULL)+
  ylab('OS(Months)')+
  
  # theme_classic()+
  theme(
    legend.position = "none",
    axis.text.x.bottom = element_blank(),  #修改坐标轴文本大小
    axis.ticks = element_blank(),   #不显示坐标轴刻度
    legend.title = element_blank(),   #不显示图例title
    # 设置背景参数
    # panel.grid.major = element_blank(),  # 去除背景网点，下面也是
    # panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#f3f3f3"),   # 设置背景颜色
    # panel.background = element_blank(),
    # 设置坐标轴与边框参数
    axis.line = element_line(size = 0.3),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
  )
p_Scatter

# p0  # 是radscore条
# p_Scatter # 散点图

# p_group
# p_baseline
# p_pathology
# p_location
# p_rna


ret <- p0/p_Scatter/p_group/p_baseline/p_pathology/p_location/p_rna +
  plot_layout(heights = c(1,
                          12,
                          1,1,5,5,1.5)
              )
ret




ggsave('复合热图-test.pdf', ret, width = 10.8, height = 7,dpi = 300)#保存文件





