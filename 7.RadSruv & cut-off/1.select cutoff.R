rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)

# █████████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1.数据准备
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
df <- read.csv('20_RiskScore_radiomics1402_pretreated_OS.csv', header=T,row.names = 1)

# 筛出训练集数据
TCGA_LGG <- read.csv('../分组/tcga_lgg.txt', header=F)
TCGA_GBM <- read.csv('../分组/tcga_gbm.txt', header=F) 
RBRT_LGG <- read.csv('../分组/rbrt_lgg.txt', header=F)
RBRT_GBM <- read.csv('../分组/rbrt_gbm.txt', header=F)
JSPH_LGG <- read.csv('../分组/jsph_lgg.txt', header=F) 
JSPH_GBM <- read.csv('../分组/jsph_gbm.txt', header=F)

list_all = rbind(TCGA_LGG,TCGA_GBM,RBRT_LGG,RBRT_GBM,JSPH_LGG,JSPH_GBM)$V1
list_lgg = rbind(TCGA_LGG,RBRT_LGG,JSPH_LGG)$V1
list_gbm = rbind(TCGA_GBM,RBRT_GBM,JSPH_GBM)$V1

df_all = df[list_all,]
df_lgg = df[list_lgg,]
df_gbm = df[list_gbm,]

rm(df,
   TCGA_LGG,TCGA_GBM,JSPH_LGG,JSPH_GBM,RBRT_LGG,RBRT_GBM,
   list_all,list_lgg,list_gbm)


# 计算中位数及截断值的函数
find_best_cutoff <- function(df_demo){
  # df_demo = df_train_all
  # 查看中位数
  RF_median = median(df_demo$RiskScore_h0)
  cat('中位数RF为:',RF_median,',')
  
  # 查看平均数
  RF_mean = mean(df_demo$RiskScore_h0)
  cat('平均数RF为:',RF_mean,',')
  
  # 取surv_cutpoint函数值做截断值？
  res.cut <- surv_cutpoint(df_demo, #数据集
                           time = "time", #生存状态
                           event = "status", #生存时间
                           variables = c("RiskScore_h0") #需要计算的数据列名
  )
  RF_cut = summary(res.cut)$cutpoint #查看数据最佳截断点及统计量

  cat('截断值RF为:',RF_cut,'\n')
  return(cbind(RF_cut,RF_median,RF_mean))
}



# 设置100次迭代，根据LGG和GBM找到最佳的cutoff值
iteration = 500
final_ret = cbind()
for (itr in 1:iteration){
  # 分割方案
  ratio = sample(c(0.3,0.4,0.5,0.6,0.7),1)
  
  sub_all<-df_all[sample(1:nrow(df_all),round(nrow(df_all)*ratio)),]
  sub_lgg<-df_lgg[sample(1:nrow(df_lgg),round(nrow(df_lgg)*ratio)),]
  sub_gbm<-df_gbm[sample(1:nrow(df_gbm),round(nrow(df_gbm)*ratio)),]
  
  cat('————————————————————————————————————————\n',
      "当前分割",ratio,'\n','lgg',nrow(sub_lgg),',gbm',nrow(sub_gbm),',lgg+gbm',nrow(sub_all),'\n')
  
  cutoff_lgg <- find_best_cutoff(sub_lgg)
  cutoff_gbm <- find_best_cutoff(sub_gbm)
  cutoff_all <- find_best_cutoff(sub_all)
  
  ret <- cbind(cutoff_lgg,cutoff_all,cutoff_gbm)
  ret <- as.data.frame(ret)
  
  # 给列命名
  names(ret) <- c('cut_off_lgg','median_lgg','mean_lgg',
                  'cut_off_all','median_all','mean_all',
                  'cut_off_gbm','median_gbm','mean_gbm')
  
  final_ret <- rbind(final_ret,ret,stringsAsFactors = F)
  
}
final_ret


# # 删掉cutoff值中特别大的，比如大于3.5(显然是有问题)
# final_ret = final_ret[-which(final_ret$cut_off_gbm > 3.5),]



# 保存上述的值
median_cut_off = as.data.frame(rbind(
  median(final_ret$cut_off_lgg),
  median(final_ret$cut_off_all),
  median(final_ret$cut_off_gbm),
  median(final_ret$median_lgg),
  median(final_ret$median_all),
  median(final_ret$median_gbm),
  median(final_ret$mean_lgg),
  median(final_ret$mean_all),
  median(final_ret$mean_gbm)
  ))
rownames(median_cut_off) <- c('cut_off_lgg','cut_off_all','cut_off_gbm',
                              'median_lgg','median_all','median_gbm',
                              'mean_lgg','mean_all','mean_gbm')
write.table(median_cut_off,'→median_and_cutoff.txt',sep = ",",quote =F,row.names=T,col.names = F)
median_cut_off


# 变长数据
dfLong <- final_ret %>%
  rownames_to_column("name1") %>%
  pivot_longer(-1,names_to = "class",values_to = "Value")
dfLong <- dfLong[,c('class','Value')]


# 新加入一列，表示是midian 还是cutoff
dfLong$cutoff_type <- dfLong$class
for (i in 1:nrow(dfLong)){
  cell = dfLong$cutoff_type[i]
  if (grepl('median',cell)){
    dfLong$cutoff_type[i] <- 'median'
  }
  if (grepl('mean',cell)){
    dfLong$cutoff_type[i] <- 'mean'
  }
  if (grepl('cut_off',cell)){
    dfLong$cutoff_type[i] <- 'cut_off'
  }
}


# 调整顺序
dfLong$class <- factor(dfLong$class,
                       levels = c('cut_off_lgg','cut_off_all','cut_off_gbm',
                                  'median_lgg','median_all','median_gbm',
                                  'mean_lgg','mean_all','mean_gbm')
                       # levels = rownames(df_corr)
)
# 再变成数据框
dfLong <- as.data.frame(dfLong)


# ?ggboxplot
plot = ggboxplot(dfLong, x = "class", y = "Value",
          color = "class",
          alpha = 0.1,
          fill="class",
          size=0.8,   #箱型图箱子边线的粗细
          add.params =list(size = 0.5, jitter = 0.25),

          palette = "jco",    # lancet  jco
          add = "jitter",
          # outlier.shape=NA, #不显示outlier
          facet.by = "cutoff_type")+
  coord_cartesian(ylim = c(0.5, 4.2))+
  rotate_x_text(-45)

plot


ggsave("1.cut-off值情况.pdf", plot, width = 15, height = 7,dpi = 300)#保存文件


