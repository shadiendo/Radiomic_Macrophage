rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(survival)
library(survminer)
library(ggplot2)


# █████████████████████████████████████████████████████████████████████████████
# 打开数据
# █████████████████████████████████████████████████████████████████████████████
estimate_df = read.csv('→→IHC和影像组学数据交集.csv', header=T,row.names = 1)
colnames(estimate_df)


surv_TTP<<-survfit(Surv(time, status) ~ WhichRiskBox_h0, data =estimate_df)

# 设置图例的文字内容
num_LowRisk = 0
for (i in estimate_df$WhichRiskBox_h0){if (i == 'low'){num_LowRisk=num_LowRisk+1}}
num_HighRisk = nrow(estimate_df)-num_LowRisk

legend_LowRisk = paste0('Low Risk (n=',num_LowRisk,')')
legend_HighRisk = paste0('High Risk (n=',num_HighRisk,')')

ggg <<- ggsurvplot(surv_TTP,data = estimate_df,pval = TRUE,
                   conf.int = FALSE,
                   risk.table = FALSE,
                   font.main=18,
                   xlab = "Overall Survival (Months)",  #x轴的label
                   ylab = "Survival Proportion",#y轴的label
                   
                   # 主线条
                   linetype = "strata",    # 根据分组修改线的类型
                   palette = c('#2b9ddf','#f0b401'),  # 颜色
                   # palette = 'aaas',   # Okabe-Ito
                   
                   # 坐标轴
                   surv.scale = 'percent',
                   surv.median.line = "hv",    # 把中位生存期显示
                   break.time.by = 12,  # x轴间隔，每24个月显示一格
                   legend.title = "",   # 图例不显示前面的
                   legend = c(0.80,0.9), # 指定图例位置
                   legend.labs = c(legend_LowRisk, legend_HighRisk),  # 指定图例分组标签
                   xlim = c(0, 48),
                   # 
                   # ggtheme=custom_theme(),
)
ggg


splots <- list()
splots[[1]] <- ggg

# Arrange multiple ggsurvplots and print the output
res = arrange_ggsurvplots(splots, print = FALSE,
                    ncol = 1, nrow = 1, risk.table.height = 1)

ggsave('CD138.163_survplots.pdf', res, width =6, height = 6, dpi = 300)




