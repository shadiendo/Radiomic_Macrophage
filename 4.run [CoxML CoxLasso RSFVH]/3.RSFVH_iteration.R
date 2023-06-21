rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(survival)
library(randomForestSRC)
library(Matrix)
library(glmnet)
library(survminer)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1.数据准备
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
TCGA_LGG <- read.csv('../分组/tcga_lgg.txt', header=F)
TCGA_GBM <- read.csv('../分组/tcga_gbm.txt', header=F) 
RBRT_LGG <- read.csv('../分组/rbrt_lgg.txt', header=F)
RBRT_GBM <- read.csv('../分组/rbrt_gbm.txt', header=F)
JSPH_LGG <- read.csv('../分组/jsph_lgg.txt', header=F) 
JSPH_GBM <- read.csv('../分组/jsph_gbm.txt', header=F)

LGG_list = rbind(TCGA_LGG,JSPH_LGG,RBRT_LGG)
GBM_list = rbind(TCGA_GBM,JSPH_GBM,RBRT_GBM) 

# 清除变量
rm(TCGA_LGG,TCGA_GBM,JSPH_LGG,JSPH_GBM,RBRT_LGG,RBRT_GBM) 

file_name <- '→radiomics1402_pretreated_OS.csv'
df_orig <- read.csv(file_name, header=T, row.names=1)   # 打开文件  , row.names=1
df_orig = df_orig[row.names(df_orig) %in% rbind(GBM_list,LGG_list)$V1,]  # 选择被试
df_orig <- as.data.frame(lapply(df_orig,as.numeric),row.names=rownames(df_orig))    # 将数据类型全部转成double

# 创建要处理的特征名,删除不需要的,保留要处理的特征 X ,交给下面的生存分析
X <- names(df_orig)     
for(item in c("time","status")){X=X[-which(X==item)]};rm(item);


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1.主函数运行
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
mianFunction <- function(rto,itr){
  # ratio = 0.3  # 分割比例
  # iteration = 1   # 迭代次数
  ratio = rto  # 分割比例
  iteration = itr   # 迭代次数
  
  # 这里预留一步，建立所有特征对应的计数表，用于做后续的迭代
  statistics <- data.frame(featureName = X,parm = 0,parm_lgg = 0,parm_gbm = 0,parm_conjugate=0)
  # 文件保存的名称
  saveFile <- paste0(iteration,'iRSFVH_',rto,'.csv')
  
  # 迭代1000次，进行统计
  for (kkk in 1:iteration){
    sub<-sample(1:nrow(df_orig),round(nrow(df_orig)*ratio))  # 分割方案
    df_train <- df_orig[sub,] #训练集
    df_test <- df_orig[-sub,] #测试集
    cat("当前的分割为",ratio,',训练集',nrow(df_train),'行,测试集',nrow(df_test),'行\n')
    
    # 分出测试集中的高级别和低级别
    df_test_GBM = na.omit(df_test[GBM_list$V1,])
    df_test_LGG = na.omit(df_test[LGG_list$V1,])
    cat("————→ 测试集中,高级别",nrow(df_test_GBM),'行,低级别',nrow(df_test_LGG),'行\n')
    
    rm(sub)
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 2.1【训练集上】构建单因素COX回归模型
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    univ_formulas <- sapply(X,function(x) as.formula(paste('Surv(time, status)~', x)))
    univ_models <- lapply( univ_formulas, function(x){coxph(x, data = df_train)})
    
    # 提取数据，并制作数据表格 
    univ_results <- lapply(univ_models,
                           function(x){
                             x <- summary(x)
                             p.value<-signif(x$wald["pvalue"], digits=2) # signif四舍五入，digits保留数字位数为2
                             wald.test<-signif(x$wald["test"], digits=2)
                             beta<-signif(x$coef[1], digits=2);    # coeficient beta
                             HR <-signif(x$coef[2], digits=2);    # exp(beta)
                             HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                             HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                             HR <- paste0(HR, " (", 
                                          HR.confint.lower, "-", HR.confint.upper, ")")
                             res<-c(beta, HR, wald.test, p.value)
                             names(res)<-c("coef", "HR (95% CI)", "wald.test", 
                                           "p.value")
                             return(res)
                           })
    univ_resTable <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
    
    # 清洗
    rm(univ_formulas,univ_models,univ_results)
    
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 2.2【训练集】筛选掉 p>0.05 的特征
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Na_pValue = is.na(univ_resTable$p.value)   # p.value是否为空的结果
    Na_Coef = is.na(univ_resTable$coef)   # coef是否为空的结果
    SumOfNa_1 = sum(Na_pValue)
    SumOfNa_2 = sum(Na_Coef)
    if (SumOfNa_1 != 0 | SumOfNa_2!=0){
      cat("当前单因素回归结果中p.value有",SumOfNa,',个缺失值，请检查？\n')
      # 排除掉所有单因素回归中是空的结果（为什么会空？可能是没有做好预处理，比如一列全是1）
      univ_resTable = univ_resTable[-which(Na_pValue == TRUE | Na_Coef==TRUE),]  
    }
    
    # 修改 p.value,coef 列的数据类型
    for (item in c('p.value','coef')){univ_resTable[,item] = as.numeric(univ_resTable[,item])};rm(item); 
    # sapply(univ_resTable, typeof)
    
    # 筛选p值
    univ_resTable <- univ_resTable[-which(univ_resTable$p.value > 0.05),]
    
    # 如果没有有意义的行数，直接退出当前循环 
    if (nrow(univ_resTable) ==0 ){next}
    
    
    # 根据筛选，提取出原列表中目标列
    df_screened<- df_train[ ,colnames(df_train) %in% row.names(univ_resTable)]
  
    # 合并一下OS信息
    df_screened <- cbind(df_train[c('time','status')],
                         df_train[ ,colnames(df_train) %in% colnames(df_screened)])
    
    # 清洗
    rm(Na_pValue,Na_Coef,SumOfNa_1,SumOfNa_2)
    
    
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # https://zhuanlan.zhihu.com/p/215883718
    # 2.3【训练集】构建RFSVH模型,最小深度法特征选择
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 模型
    # nodesize_best = tune.nodesize(Surv(time, status) ~ ., df_screened)
    # nodesize_best = nodesize_best$nsize.opt
    nodesize_best = 3
    ret.obj <- rfsrc(Surv(time, status) ~ ., data=df_screened, 
                     nodesize = nodesize_best, 
                     importance = TRUE)
    # # 绘图
    # plot(ret.obj)
    # plot(get.tree(ret.obj, 3))
    
    # different levels of conservativeness,  high medium  low
    ret.vb <- var.select(object = ret.obj, conservative = "high")
    
    
    allvars_condition <- ret.vb$varselect   # 每个特征的评分表,分别为depth和vimp指标
    meaningful_fn <- ret.vb$topvars
    
    
    # # 尝试删除vimp为负，且vimp<0.01的最好特征
    # meaningful_fn = NA
    # for (i in topvars){
    #   if (allvars_condition[i,2] > 0.005){
    #     meaningful_fn = rbind(meaningful_fn,i)
    #   }
    # }
    # meaningful_fn = as.data.frame(meaningful_fn)$V1
    # meaningful_fn = meaningful_fn[2:length(meaningful_fn)]
    # meaningful_fn
    
    
    if (length(meaningful_fn)==0){next}
    
  
    # 清理
    rm(nodesize_best,ret.obj,ret.vb,allvars_condition)
  
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 定义两个函数
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    predictTheOS <- function(df_demo){
      # 根据上述筛选出来的特征
      df_demo$RF <- 0
      df_demo$RF_whichBox <- 0
      
      for (eachPat in 1:nrow(df_demo)){
        # 对每一行计算RF值 RF值=∑featureValue*ceof
        RF_for_this_one = 0     
        
        for (col in meaningful_fn){
          cell_value <- df_demo[eachPat,col]
          cell_value[is.na(cell_value)] <- 0   # 如果值为NA，化为0
          
          eachRF =  cell_value*univ_resTable[col,'coef']   
          # cat('(',cell_value,'*',res[col,'coef'],')+')   # 打印出公式
          
          RF_for_this_one = RF_for_this_one+eachRF}  # 加入该病人的计算公式
        
        # 将该病人RF得分加入新数据框data_new中
        df_demo[eachPat,'RF'] = RF_for_this_one
      }
      
      # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
      # 1 计算每个被试的RF评分，求得平均RF
      # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
      # 取中位数做截断值？
      RF_median = median(df_demo$RF)
      cat('中位数RF为:',RF_median,'\n')
      
      # 取surv_cutpoint函数值做截断值？
      res.cut <- surv_cutpoint(df_demo, #数据集
                               time = "time", #生存状态
                               event = "status", #生存时间
                               variables = c("RF") #需要计算的数据列名
      )
      RF_cut = summary(res.cut)$cutpoint #查看数据最佳截断点及统计量
      
      RF_cut = RF_cut
      cat('截断值RF为:',RF_cut,'\n')
      
      # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
      # 2 根据平均RF，拆分两组，低于RF_cut标记为0，高于则标记为1
      # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
      for (eachPat in 1:nrow(df_demo)){
        cell_value <- df_demo[eachPat,'RF']
        if (cell_value<RF_cut){
          df_demo[eachPat,'RF_whichBox'] <- 0
        }else{
          df_demo[eachPat,'RF_whichBox'] <- 1
        }
      }
      
      # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
      # 3 Kaplan-Meier生存分析,根据log rank检的P值,确定这组特征是否保留
      # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
      surv_diff <- survdiff(Surv(time,status) ~ RF_whichBox,data = df_demo)
      pvalue = surv_diff$pvalue
      # cat('KM生存分析 P值:',pvalue,'(log rank检验)\n')
      
      return(pvalue)
      
      # 清洗
      rm(col,eachPat,RF_for_this_one,cell_value,eachRF,
         RF_mean,
         surv_diff,pvalue
      )
    }
    
    writeOrNot <- function(pvalue_demo,colname_demo){
      if (pvalue_demo<0.05){
        for (fname in meaningful_fn){
          # 定位到单元格[fname,parm]
          cnt <- statistics[which(statistics$featureName == fname),colname_demo]+1   
          # 写入计数,【注意!!!】这里因为是函数内，要修改函数外的东西，得使用  <<- 
          statistics[which(statistics$featureName == fname),colname_demo] <<- cnt}
      }
    }
    
    
    
    # 训练集：如果这轮特征连自己都满足不了，那么下一轮
    pvalue.train = predictTheOS(df_train)
    cat('pvalue.train:',pvalue.train,'\n')
    if (pvalue.train>=0.05){next}
    
    # 测试集
    pvalue.test = predictTheOS(df_test)
    writeOrNot(pvalue.test,'parm')
    cat('pvalue.test',pvalue.test,'\n')
    
    # 测试集-lgg
    pvalue.test.lgg = predictTheOS(df_test_LGG)
    writeOrNot(pvalue.test.lgg,'parm_lgg')
    cat('pvalue.test.lgg',pvalue.test.lgg,'\n')
    
    # 测试集-gbm
    pvalue.test.gbm = predictTheOS(df_test_GBM)
    writeOrNot(pvalue.test.gbm,'parm_gbm')
    cat('pvalue.test.gbm',pvalue.test.gbm,'\n')
    
    
    # 这些特征是否在 LGG&GBM LGG GBM 都有区分意义
    if (pvalue.test<0.05 && pvalue.test.lgg<0.05 && pvalue.test.gbm<0.05){
      cat('本轮特征具有同时区分高低级别，高级别和低级别的性能')
      for (fname in meaningful_fn){
        cnt <- statistics[which(statistics$featureName == fname),'parm_conjugate']+1
        statistics[which(statistics$featureName == fname),'parm_conjugate'] <- cnt
      }
    }
    
  
    
    cat('\n┏━━━━━━━━━━━━━━━━╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱━━━━━━━━━━━━━━━━┓',  
        '\n┃                iteration:',kkk,'finished.                 ',
        '\n┗━━━━━━━━━━━━━━━━╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱━━━━━━━━━━━━━━━━┛\n\n')
    
  }
  
  # 排序并改名
  # statistics <- statistics[order(-statistics[,2]),]   
  names(statistics) <- c('featureName',
                         paste0('parm',rto*10,(1-rto)*10),
                         paste0('parm',rto*10,(1-rto)*10,'_lgg'),
                         paste0('parm',rto*10,(1-rto)*10,'_gbm'),
                         paste0('parm',rto*10,(1-rto)*10,'_conjugate'))
  
  # 保存文件
  write.csv(statistics,saveFile)
}

mianFunction(0.3,1000)
mianFunction(0.4,1000)

