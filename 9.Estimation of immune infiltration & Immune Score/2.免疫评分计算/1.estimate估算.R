# library(utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge, dependencies=TRUE)

rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(estimate)


# # ——————————————————————————————————————————————————
# # 打开示例数据（择一）
# # ——————————————————————————————————————————————————
# OvarianCancerExpr <- system.file("extdata", "sample_input.txt", package="estimate")
# sample <- read.csv(OvarianCancerExpr, header=T,sep="\t")

# ——————————————————————————————————————————————————
# 打开自己的数据（择一）
# ——————————————————————————————————————————————————
# 打开自己的数据转成txt
open_file = 'TCGA/TCGA_LGG_fpkm.csv'
df_orig <- read.csv(open_file, header=T)   # 打开文件 row.names=1
write.table(df_orig,file="exp.txt",row.names=FALSE,sep="\t",quote = F)

OvarianCancerExpr <- "exp.txt"
sample <- read.csv(OvarianCancerExpr, header=T,sep="\t")


# ——————————————————————————————————————————————————
# 运行
# ——————————————————————————————————————————————————
# https://www.jianshu.com/p/ec5307256ca5
# http://t.csdn.cn/FuuIA
# 中间步骤
filterCommonGenes(input.f=OvarianCancerExpr, 
                  output.f="OV_10412genes.gct", 
                  id="GeneSymbol")

#进行estimate分析
estimateScore(input.ds="OV_10412genes.gct",
              output.ds="OV_estimate_score.gct", 
              platform="affymetrix")
# platform = c("affymetrix", "agilent", "illumina"))


# # 手动画图
# plotPurity(scores="OV_estimate_score.gct", 
#            # samples="s519",   # 如果注释掉这条参数的话，就是把所有的图一次性画出来
#            platform="affymetrix",
#            output.dir="estimated_purity_plots")


# 将gct文件整理成csv
estimate_score <- read.table("OV_estimate_score.gct", skip = 2, header = TRUE)
rownames(estimate_score) <- estimate_score[,1]
estimate_score <- estimate_score[,3:ncol(estimate_score)]


saveFile <- paste0(getwd(),'/',substr(open_file,1,nchar(open_file)-4),'_estimate_ret.csv') 
write.csv(estimate_score,saveFile)
cat('已保存至',saveFile)



# 删除中间文件
file.remove("exp.txt","OV_10412genes.gct","OV_estimate_score.gct")


# 
# # ——————————————————————————————————————————————————
# # 绘图?
# # ——————————————————————————————————————————————————
# #生成plot
# scores=read.table("OV_estimate_score.gct",skip = 2,header = T)
# rownames(scores)=scores[,1]
# scores=t(scores[,3:ncol(scores)])
# View(scores)
# scores<- as.data.frame(scores)
# scores$SampleID <- rownames(scores)
# 
# save(scores,file = 'BRCA_estimate_score.rdata') 
# # m_score <- get(load('BRCA_estimate_score.rdata'))
# class(scores)
# 
# # median(scores[,'ImmuneScore'])
# 
# ## 画单条累计分布曲线
# library(ggplot2)
# p1<-ggplot(scores,aes(x=ImmuneScore)) +
#   stat_ecdf(color = "green") +
#   labs(y="accumulative propotion")
# 
# p2<-ggplot(scores,aes(x=StromalScore)) +
#   stat_ecdf(color = "red") +
#   labs(y="accumulative propotion")
# 
# 
# ## 画多条累计分布曲线
# # 要先把数据框转化成长数据格式
# library(reshape2)
# scores2 <- scores[,-4] # 去掉不需要的列
# dim(scores2)
# long = melt(scores2, id=c("SampleID"),
#             variable.name= 'Class', value.name = 'Value')
# # 创建新的变量，含有分组信息
# long$group <- rep(c(rep("g1",6),rep("g2",4)),3)
# library(dplyr)
# long<- mutate(long,group_score = paste(Class,group,sep="_"))
# 
# ## 长格式转成宽格式数据框
# #wide <- dcast(long,SampleID~Class,value.var='Value')
# 
# #作多条累计分布曲线，以不同颜色区分。
# ggplot(long,aes(x=Value,color =Class,linetype=group)) +
#   stat_ecdf(size =1) + # 线粗细
#   labs(x="Score",
#        y="Accumulative propotion",
#        title="Accumulative Plotting") +
#   theme(legend.position = c(0.15,0.7), #图例位置 ("none", "left", "right", "bottom", "top", or two-element numeric vector)
#         legend.background =element_rect(fill = "white", colour = "grey50"),#图例背景
#         panel.background=element_rect(fill = "grey95", # grey90
#                                       colour = "black",
#                                       size = 1), #画布背景颜色
#         plot.title = element_text(hjust=0.5,size=16,vjust=0.5), #标题位置
#         legend.text=element_text(size=10,colour='black'), #图例文字
#         axis.text=element_text(size=8,colour="black"), #坐标轴文字
#         axis.title.y = element_text(size = rel(1.3), angle = 90),#y坐标轴名称文字
#         axis.title.x = element_text(size = rel(1.3)),#x坐标轴名称文字
#   )
