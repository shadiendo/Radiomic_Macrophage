rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(ggplot2)
library(ggpubr)
library(patchwork)


# 打开数据集
open_file <- '→长数据_供绘散点图.csv'
df_orig <- read.csv(open_file, header=T)  # ,row.names=F
df_orig$riskScore_log = log2(df_orig$RiskScore_h0)  # 新增一列 riskScore_log

# 要不要删除所有用不到的
MacrophageList = cbind('StromalScore','ImmuneScore','ESTIMATEScore')
df_orig = df_orig[df_orig$class %in% MacrophageList,]



# 筛出训练集、测试集数据
TCGA_LGG <- read.csv('../../分组/tcga_lgg.txt', header=F)$V1
TCGA_GBM <- read.csv('../../分组/tcga_gbm.txt', header=F)$V1
TCGA <- c(TCGA_LGG,TCGA_GBM)

RBRT_LGG <- read.csv('../../分组/rbrt_lgg.txt', header=F)$V1
RBRT_GBM <- read.csv('../../分组/rbrt_gbm.txt', header=F)$V1
RBRT <- c(RBRT_LGG,RBRT_GBM)

JSPH_LGG <- read.csv('../../分组/jsph_lgg.txt', header=F)$V1
JSPH_GBM <- read.csv('../../分组/jsph_gbm.txt', header=F)$V1
test_JSPH_LGG <- read.csv('../../分组/TEST_jsph_lgg.txt', header=F)$V1
test_JSPH_GBM <- read.csv('../../分组/TEST_jsph_gbm.txt', header=F)$V1
JSPH <-c(JSPH_LGG,JSPH_GBM,test_JSPH_LGG,test_JSPH_GBM)

LGG <- c(TCGA_LGG,RBRT_LGG)
GBM <- c(TCGA_GBM,RBRT_GBM,JSPH_GBM)

ALL <- c(TCGA,RBRT)
ALL_LGG <- c(TCGA_LGG,RBRT_LGG,JSPH_LGG,test_JSPH_LGG)
ALL_GBM <- c(JSPH_GBM,test_JSPH_GBM,RBRT_GBM,TCGA_GBM)

# ████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
showThePlot <- function(group,title_name) {
  # group = TCGA
  # title_name = 'TCGA'
  
  df_demo = df_orig[which(df_orig$X %in% group),]
  
  gg <- ggscatter(df_demo, x = "enrichment_value", y = "riskScore_log",
                  color='class',
                  # size = "Correlation",
                  # ———————————————————————————————————————————
                  add = "reg.line",
                  conf.int = TRUE,
                  # add.params = list(fill = "lightgray"),
                  # ———————————————————————————————————————————
                  palette = "lancet",     # 配色方案杂志 jama,nature,lancet,jco
                  # ———————————————————————————————————————————
                  # ellipse = TRUE,  # 画椭圆
                  # mean.point = TRUE,  # 中值
                  # star.plot = TRUE   #生成星图
                  # ———————————————————————————————————————————
  )+
    # scale_y_continuous(limits =c(-4, 3))+
    # scale_x_continuous(limits =c(0, 0.6))+
    
    stat_cor(aes(color = class),
             show.legend = FALSE,
             method = "pearson",
             size=5,
             r.accuracy = 0.001,
             p.digits = 3,
             # # label.x.npc="left",
             label.x = 0.4,
             label.y.npc = 0.2,
             # # label.y = 0
    )
  gg
  
  gg = ggpar(gg,
             title = title_name,
             # title = paste0("Correlation Scatter Diagram (",title_name,')'),
             # subtitle = "Relationship between macrophages and RiskScore",
             font.title = c(20, "bold","black"),
             font.subtitle = c(12, "italic","black"),
             
             xlab = "Macrophage_enrichment",
             ylab = "log2(RiskScore)",
             font.x = c(12, "black"),
             font.y = c(12, "black"),
             
             legend.title = "",
             # legend = "right"
  )
  gg
  
}

plot_TCGA_LGG = showThePlot(TCGA_LGG,'TCGA_LGG')
plot_TCGA_GBM = showThePlot(TCGA_GBM,'TCGA_GBM')
plot_TCGA = showThePlot(TCGA,'TCGA')

plot_RBRT_LGG = showThePlot(RBRT_LGG,'RBRT_LGG')
plot_RBRT_GBM = showThePlot(RBRT_GBM,'RBRT_GBM')
plot_RBRT = showThePlot(RBRT,'RBRT')

plot_JSPH = showThePlot(JSPH,'JSPH')

plot_LGG = showThePlot(LGG,'LGG')
plot_GBM = showThePlot(GBM,'GBM')

plot_ALL = showThePlot(ALL,'ALL')
plot_ALL_LGG = showThePlot(ALL_LGG,'ALL_LGG')
plot_ALL_GBM = showThePlot(ALL_GBM,'ALL_GBM')


plot_ret <- (plot_TCGA_LGG | plot_TCGA_GBM | plot_TCGA | plot_ALL)/
  (plot_RBRT_LGG | plot_RBRT_GBM | plot_RBRT | plot_ALL_LGG)/
  (plot_LGG | plot_GBM | plot_JSPH| plot_ALL_GBM)
ggsave('→散点图结果汇总.pdf', plot_ret, width = 28, height = 21,dpi = 300)#保存文件






