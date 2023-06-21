rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(ggplot2)
library(ggpubr)
library(patchwork)


# 打开数据集
open_file <- 'demo_data.csv'
df_orig <- read.csv(open_file, header=T,encoding="UTF-8")  # ,row.names=F
colnames(df_orig)


showThePlot <- function(group) {
  ggboxplot(df_orig, x = "WhichRiskBox_h0", y = group,
            color = "WhichRiskBox_h0", 
            palette = "jco",
            add = "jitter",
            title = group,
  )+
    stat_compare_means(method = "t.test")   # anova  t.test wilcox.test
}


CD20.CD68 = showThePlot('CD20.CD68')
CD138.163 = showThePlot('CD138.163')
# IgG.PTEN = showThePlot('IgG.PTEN')


plot_ret <- (CD20.CD68 | CD138.163)
plot_ret


ggsave('→箱线图结果汇总.pdf', plot_ret, width = 8, height = 6,dpi = 300)#保存文件


# t.test(CD20.CD68_mean_value~risklabel_new, var.equal=T,data=df_orig)



