rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(survival)
library(survminer)
library(ggplot2)
# remove.packages('purrr')
# install.packages('purrr')


riskScore_manually <- function(df_file_name,multi_resTable_file_name) {
  
  # 读取影像组学OS文件
  df_demo = read.csv(df_file_name, header=T,row.names = 1)
  
  # 读取多因素风险文件
  multi_resTable <- read.csv(multi_resTable_file_name, header=T,row.names = 1)
  feature_mutiCox = rownames(multi_resTable)
  
  # 设置保存的名称
  save_name = paste0('→→→RiskScore_',df_file_name)
  
  
  # 计算未经修正的riskScore用于riskScore_h0的基础
  df_demo$RiskScore <- 0
  for (eachPat in 1:nrow(df_demo)){
    # 对每一行计算RF值 RF值=∑featureValue*ceof
    RiskScore_value = 0
    
    for (col in feature_mutiCox){
      cell_value <- df_demo[eachPat,col]
      cell_value[is.na(cell_value)] <- 0   # 如果值为NA，化为0
      
      eachRF =  cell_value*multi_resTable[col,'coef']
      # cat('(',cell_value,'*',res[col,'coef'],')+')   # 打印出公式
      
      RiskScore_value = RiskScore_value+eachRF}  # 加入该病人的计算公式
    
    # 将该病人RF得分加入新数据框data_new中
    df_demo[eachPat,'RiskScore'] = RiskScore_value
  }
  
  # -----------------------------------------------------------------------
  df_demo$RiskScore_h0 <- 0
  # 创建h0修正后的RiskScore
  for (eachPat in 1:nrow(df_demo)){
    h0_corrected <- exp(h0+df_demo[eachPat,'RiskScore'])
    df_demo[eachPat,'RiskScore_h0'] <- h0_corrected
  }
  
  # -----------------------------------------------------------------------  
  df_demo <- cbind(df_demo[c('RiskScore_h0')],
                   df_demo[ ,colnames(df_demo) %in% feature_mutiCox])
  

  # 保存RiskScore的结果
  write.csv(df_demo,save_name)
  cat('【RF计算结果保存成功!!!】\n')
}



h0= -206.4363
cat('当前训练数据的h0为',h0)

riskScore_manually(
  'mass.csv',
  '→→mutiCox_coef_num20.csv')





