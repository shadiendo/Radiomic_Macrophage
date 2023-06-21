rm(list=ls())

library(this.path)
setwd(this.path::this.dir())
cat('当前脚本执行路径为',getwd())

source("source/cibersort.R")  # 引用

# 打开数据集
sig_matrix <- "source/sig_matrix.txt"   # 注释文件名
mixture_file <- '→JSPH_exp_B.csv'

# # 这里其实可以不用打开，因为cibersort源码识别文件名
# df_orig <- read.csv(mixture_file, header=T,row.names=1)   # 打开文件 row.names=1

# -----------------------------------
# 【运行CIBERSORT代码】
# 注：
# CIBERSORT源码打开的矩阵文件是txt： read.table(mixture_file, header=T, sep="\t", check.names=F)
# 如果你要打开csv,你可以手动修改成： read.csv(mixture_file, header=T)

res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)  #用cibersort算法计算，perm=100模拟计算次数,决定p值
res_cibersort = as.data.frame(res_cibersort)

save_name <- paste0('★',substr(mixture_file,1,nchar(mixture_file)-4),'_cibersort.csv')
cat('保存文件名为',save_name)
write.csv(res_cibersort,save_name)



# #############################################
#可视化展示
# #############################################
rm(list=ls())

# open_filename = save_name
open_filename = '★Rembrandt_expression_trans_cibersort.csv'

df_orig <- read.csv(open_filename, header=T,row.names=1)   # 打开文件 row.names=1
res_cibersort <- df_orig[,1:22]   #取前22列为细胞丰度数据
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #去除丰度全为0的细胞

#barplot图
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
a = barplot(as.matrix(t(ciber.res)),
            border = NA, # 柱子无边框
            names.arg = rep("",nrow(ciber.res)), # 无横坐标样本名
            yaxt = "n", # 先不绘制y轴
            ylab = "Relative percentage", # 修改y轴名称
            col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-0,
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.8, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   #关闭画板

#相关性热图
#install.packages('corrplot')
M <- round(cor(ciber.res),2) # 计算相关性矩阵并保留两位小数

library(corrplot)
corrplot.mixed(M,
               lower.col = "black", #左下方字体颜色为黑色
               tl.pos = "lt",  #标签出现在左侧和顶部
               number.cex = 0.5, #左下方字号为0.5
               tl.cex = 0.5) #标签字号为0.5
dev.off()   #关闭画板


#各个细胞组分的含量箱线图
library(tidyr)
library(dplyr)
library(tidyverse)

dd1 <- ciber.res %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = 2:16,     # 【这里需要适当改一下，因为有的细胞类型结果在前面被筛掉了】
               names_to = "CellType",
               values_to = "Composition")
plot.info <- dd1[]          


library(ggplot2)
library(ggpubr)
library(ggthemes)
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Cell composition",
  main = "TME Cell composition") +
  theme_base() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 1
  ))     
