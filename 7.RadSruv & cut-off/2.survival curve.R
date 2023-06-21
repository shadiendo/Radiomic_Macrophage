rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(survival)
library(survminer)
library(ggplot2)


# █████████████████████████████████████████████████████████████████████████████
# 打开数据
# █████████████████████████████████████████████████████████████████████████████
df <- read.csv('20_RiskScore_radiomics1402_pretreated_OS.csv', header=T,row.names = 1)
save_name <- '→Compand_survplots-surv.pdf'
cut_off_lgg = 0.63628810153324
cut_off_all = 1.06124738986942
cut_off_gbm = 1.7981041963824

# 筛出训练集、测试集数据
TCGA_LGG <- read.csv('../分组/tcga_lgg.txt', header=F)
TCGA_GBM <- read.csv('../分组/tcga_gbm.txt', header=F) 
RBRT_LGG <- read.csv('../分组/rbrt_lgg.txt', header=F)
RBRT_GBM <- read.csv('../分组/rbrt_gbm.txt', header=F)
JSPH_LGG <- read.csv('../分组/jsph_lgg.txt', header=F) 
JSPH_GBM <- read.csv('../分组/jsph_gbm.txt', header=F)
test_JSPH_LGG <- read.csv('../分组/TEST_jsph_lgg.txt', header=F)
test_JSPH_GBM <- read.csv('../分组/TEST_jsph_gbm.txt', header=F)

train <- rbind(TCGA_LGG,TCGA_GBM,RBRT_LGG,RBRT_GBM,JSPH_LGG,JSPH_GBM)
train_lgg <- rbind(TCGA_LGG,RBRT_LGG,JSPH_LGG)
train_gbm <- rbind(TCGA_GBM,RBRT_GBM,JSPH_GBM)

showThePlot <- function(group,thisPlotTitle,cut_off,lowrisk_color,highrisk_color) {
  calcMedianSurvivalTime <- function(data) {
    km_fit <- survfit(Surv(time, status) ~ 1, data = data)
    
    # 提取生存概率和时间点
    surv_prob <- km_fit$surv
    time_points <- km_fit$time
    
    # 计算中位生存期
    median_survival_time <- NA
    for (i in 1:length(surv_prob)) {
      if (surv_prob[i] <= 0.5) {
        median_survival_time <- time_points[i]
        break
      }
    }
    
    return(median_survival_time)
  }
  # group = RBRT_GBM
  # thisPlotTitle = 'ret'
  # cut_off=1.553446325
  # lowrisk_color = "#a40a3c"
  # highrisk_color = "#8f8787"
  
  # 先提取
  estimate_df <<- na.omit(df[group$V1,])
  
  
  # 根据平均RF，拆分两组，低于RF_mean标记为0，高于则标记为1
  estimate_df$whitchBox <- 0
  for (eachPat in 1:nrow(estimate_df)){
    cell_value <- estimate_df[eachPat,'RiskScore_h0']
    if (cell_value<cut_off){
      estimate_df[eachPat,'whitchBox'] <- 'low'
    }else{
      estimate_df[eachPat,'whitchBox'] <- 'high'
    }
  }
  cat('该组分类中cutoff值为',cut_off,'\n')
  
  
  # cbind('RiskScore','WhichRiskBox','RiskScore_h0','WhichRiskBox_h0')
  surv_TTP<<-survfit(Surv(time, status) ~ whitchBox, data =estimate_df)
  
  # 设置图例的文字内容
  num_LowRisk = 0
  for (i in estimate_df$whitchBox){if (i == 'low'){num_LowRisk=num_LowRisk+1}}
  num_HighRisk = nrow(estimate_df)-num_LowRisk
  
  legend_LowRisk = paste0('Low Risk (n=',num_LowRisk,')')
  legend_HighRisk = paste0('High Risk (n=',num_HighRisk,')')
  legend_HighRisk
  legend_LowRisk
  
  
  ggg <<- ggsurvplot(surv_TTP,data = estimate_df,pval = TRUE,conf.int = TRUE,risk.table = FALSE,
                     # surv.plot.height=0.2,
                     title = paste0(thisPlotTitle,' (Cut-off=',round(cut_off,2),')'),
                     font.main=18,
                     xlab = "Overall Survival (Months)",  #x轴的label
                     ylab = "Survival Proportion",#y轴的label
                     
                     # 主线条
                     linetype = "strata",    # 根据分组修改线的类型
                     palette = c(lowrisk_color,highrisk_color),  # 颜色
                     # palette = 'aaas',   # Okabe-Ito
                     
                     # 坐标轴
                     surv.scale = 'percent',
                     surv.median.line = "hv",    # 把中位生存期显示
                     break.time.by = 12,  # x轴间隔，每24个月显示一格
                     legend.title = "",   # 图例不显示前面的
                     legend = c(0.80,0.9), # 指定图例位置
                     legend.labs = c(legend_HighRisk, legend_LowRisk),  # 指定图例分组标签
                     xlim = c(0, 120),
                     # 
                     # ggtheme=custom_theme(),
  )
  ggg
  
  return(ggg)
}


# cut_off_lgg = 0.995083192532716
# cut_off_all = 1.82929432512945
# cut_off_gbm = 2.86786843779624


# # 运行一个测试下
# plot_demo <- showThePlot(train_lgg,
#                          'demo',cut_off_lgg,'#a40a3c','#8f8787')
# plot_demo




# 训练集
plot_TCGA <- showThePlot(rbind(TCGA_LGG,TCGA_GBM),
                         'TCGA',cut_off_all,"#a40a3c","#8f8787")
plot_TCGA_LGG <- showThePlot(TCGA_LGG,
                             'TCGA_LGG',cut_off_lgg,"#a40a3c","#8f8787")
plot_TCGA_GBM <- showThePlot(TCGA_GBM,
                             'TCGA_GBM',cut_off_gbm,"#a40a3c","#8f8787")

plot_RBRT <- showThePlot(rbind(RBRT_LGG,RBRT_GBM),
                         'RBRT',cut_off_all,'#a40a3c','#8f8787')
plot_RBRT_LGG <- showThePlot(RBRT_LGG,
                             'RBRT_LGG',cut_off_lgg,'#a40a3c','#8f8787')
plot_RBRT_GBM <- showThePlot(RBRT_GBM,
                             'RBRT_GBM',cut_off_gbm,'#a40a3c','#8f8787')

plot_JSPH <- showThePlot(rbind(JSPH_LGG,JSPH_GBM),
                         'JSPH',cut_off_all,"#a40a3c","#8f8787")
plot_JSPH_LGG <- showThePlot(JSPH_LGG,
                             'JSPH_LGG',cut_off_lgg,"#a40a3c","#8f8787")
plot_JSPH_GBM <- showThePlot(JSPH_GBM,
                             'JSPH_GBM',cut_off_gbm,"#a40a3c","#8f8787")

# 测试集
plot_JSPH_TEST <- showThePlot(rbind(test_JSPH_LGG,test_JSPH_GBM),
                              'JSPH(TEST)',cut_off_all,'#ee0000','#a88462')
plot_JSPH_LGG_TEST <- showThePlot(test_JSPH_LGG,
                                  'JSPH_LGG(TEST)',cut_off_lgg,'#ee0000','#a88462')
plot_JSPH_GBM_TEST <- showThePlot(test_JSPH_GBM,
                                  'JSPH_GBM(TEST)',cut_off_gbm,'#ee0000','#a88462')




splots <- list()
splots[[1]] <- plot_JSPH
splots[[5]] <- plot_JSPH_LGG
splots[[9]] <- plot_JSPH_GBM

splots[[2]] <- plot_TCGA
splots[[6]] <- plot_TCGA_LGG
splots[[10]] <- plot_TCGA_GBM

splots[[3]] <- plot_RBRT
splots[[7]] <- plot_RBRT_LGG
splots[[11]] <- plot_RBRT_GBM

splots[[4]] <- plot_JSPH_TEST
splots[[8]] <- plot_JSPH_LGG_TEST
splots[[12]] <- plot_JSPH_GBM_TEST



# Arrange multiple ggsurvplots and print the output
res = arrange_ggsurvplots(splots, print = FALSE,
                          ncol = 3, nrow = 4
)

ggsave(save_name, res, width =16, height = 16, dpi = 600)







