rm(list=ls())
# install.packages("ggVennDiagram")
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(ggVennDiagram)
library(ggplot2)
# 示例数据准备


have_FLU <- read.table('▲有随访.txt',head=FALSE)$V1   # 待处理特征
have_MRI <- read.table('▲有MR.txt',head=FALSE)$V1   # 待处理特征
have_IDH <- read.table('▲有IDH.txt',head=FALSE)$V1   # 待处理特征

x<- list(A=have_FLU,
          B=have_MRI,
          C=have_IDH)


# method1
plot_ret <- ggVennDiagram(x,
              category.names = c("with OS Infomation","available MRI","with IDH status"),
              set_size = 5,
              # label = "count", 
              label_percent_digit = 1, # 百分比的位数
              label_color = "black",
              label_alpha = 0,   # 标签的背景颜色
              edge_lty = "solid", 
              edge_size = 1) +
  # legend(fill = FALSE)+
  guides(fill="none")+
  # scale_fill_distiller(palette = "PuOr", direction = -1) +
  # scale_color_brewer(palette = "Set1")
  scale_fill_gradient(low='grey95',high='grey30',name = "gene count")+   # 填充色
  scale_color_manual(values=c('black','black','black'))   # 边框色
  

ggsave('→韦恩图.pdf', plot_ret, width = 5, height = 5,dpi = 300)#保存文件


