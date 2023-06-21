rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(survival)
library(survminer)
library(ggplot2)

# █████████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1.数据准备
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
df_orig <- read.csv('radiomics1402_pretreated_OS.csv', header=T,row.names = 1)   # 影像组学数据
feature_selected <- read.table('→→muticoxMiddleRet_num20.txt',head=FALSE)$V1   # 待处理特征

# 提取目标列
df <- cbind(df_orig[c('time','status')],
                     df_orig[ ,colnames(df_orig) %in% feature_selected])
df <- as.data.frame(lapply(df,as.numeric)) 
rownames(df) <- rownames(df_orig)

rm(df_orig,feature_selected)


# 筛出训练集、测试集数据
TCGA_LGG <- read.csv('../分组/tcga_lgg.txt', header=F)
TCGA_GBM <- read.csv('../分组/tcga_gbm.txt', header=F) 
RBRT_LGG <- read.csv('../分组/rbrt_lgg.txt', header=F)
RBRT_GBM <- read.csv('../分组/rbrt_gbm.txt', header=F)
JSPH_LGG <- read.csv('../分组/jsph_lgg.txt', header=F) 
JSPH_GBM <- read.csv('../分组/jsph_gbm.txt', header=F)

train_list = rbind(TCGA_LGG,TCGA_GBM,RBRT_LGG,RBRT_GBM,JSPH_LGG,JSPH_GBM) 
df_train = na.omit(df[train_list$V1,])



# █████████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 要不要先去重？要
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# █████████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2.构建模型 - 通过训练集构建多因素COX
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
multi_Cox <- coxph(Surv(time, status) ~ ., data = df_train)
step_Cox <- step(multi_Cox ,direction = "backward")     # 参数还有forward和both
cox_summary <- summary(step_Cox)

multi_resTable = as.data.frame(cox_summary$coefficients)  # 多因素回归系数表
feature_final = rownames(multi_resTable)
cat('multi_Cox筛得',length(feature_final),'个')


# 保存特征
save_ff <- paste0('→→muticoxMiddleRet_num',length(feature_final),'.txt')
write.table(feature_final,save_ff, row.names =F, col.names =F,quote =F)
cat('暂时保存多因素回归特征')

# 保存多因素回归结果
save_mt <- paste0('→→mutiCox_coef_num',length(feature_final),'.csv')
write.csv(multi_resTable,save_mt)
cat('保存多因素回归结果csv')

rm(multi_Cox,cox_summary,feature_final,
   save_ff,save_mt)

# █████████████████████████████████████████████████████████████████████████████
# 3.计算RiskScore
# █████████████████████████████████████████████████████████████████████████████
caculate_h0 <- function(multi_resTable) {
  # https://www.jingege.wang/2022/09/26/
  # 为何predict函数计算的riskscore不等于基因的表达量与其系数
  # 
  # 确定h0
  # h(t,X) = h0(t) * exp(β1X1 + β2X2 + ... + βnXn)
  # ↓
  # lnh0(t) = lnh(t,X) - (β1X1 + β2X2 + ... + βnXn)
  # ↓
  # 即先计算predict的值，再减去手算的，为基准风险值
  
  df_demo = df_train
  feature_mutiCox = rownames(multi_resTable)

  df.demo = df_demo[1,]
  # 取第一行，算predict的值
  riskScore.predict.demo =predict(step_Cox ,type = "risk",newdata = df.demo)
  riskScore.predict.demo
  
  # 手算一个
  # multi_resTable
  riskScore.manual.demo = 0  
  for (ftr in rownames(multi_resTable)){
    # 元素值
    ftr_value <- df.demo[,ftr]
    # 元素值*HR值
    each_value =  ftr_value * multi_resTable[ftr,'coef']   
    # 打印出公式
    # cat('(',each_value,'*',multi_resTable[ftr,'coef'],')+')
    # 加入该病人的计算公式
    riskScore.manual.demo = riskScore.manual.demo+each_value
  }  
  riskScore.manual.demo
  
  # 计算h(0)
  h0 = log(riskScore.predict.demo)-riskScore.manual.demo
  h0
  return(h0)
}
h0=caculate_h0(multi_resTable)
cat('当前训练数据的h0为',h0)

# █████████████████████████████████████████████████████████████████████████████
# 3.计算所有数据 或 特定数据的RiskScore
# █████████████████████████████████████████████████████████████████████████████
# 多因素的公式为:   Riskscore = h0(t)*exp(β1X1 + β2X2 + ... + βnXn)
# R中表达为:        exp(h0+riskScore.manual.demo)

  
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
  # 查看一下中值RF
  RiskScore_median_h0 = median(df_demo$RiskScore_h0)
  cat('中值为',RiskScore_median_h0,'\n')
  
  
  # 取surv_cutpoint函数值做截断值？
  res.cut_h0 <- surv_cutpoint(df_demo, #数据集
                              time = "time", #生存状态
                              event = "status", #生存时间
                              variables = c("RiskScore_h0") #需要计算的数据列名
  )
  cut_off_best = summary(res.cut_h0)$cutpoint #查看数据最佳截断点及统计量
  cat('最佳截断值为:',cut_off_best,'\n')
  
  # -----------------------------------------------------------------------  
  df_demo <- cbind(df_demo[c('status','time','RiskScore_h0')],
                   df_demo[ ,colnames(df_demo) %in% feature_mutiCox])
  

  # 保存RiskScore的结果
  write.csv(df_demo,save_name)
  cat('【RF计算结果保存成功!!!】\n')
}



riskScore_manually(
  'radiomics1402_pretreated_OS.csv',
  '→→mutiCox_coef_num19.csv')





