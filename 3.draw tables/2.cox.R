rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(survival)
library(dplyr)
library(tableone)
library(car)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
df_orig <- read.csv('GeneralView - for cox.csv', header=T,row.names = 1,na.strings = c("", "NA"))
# 筛选行  c('JSPH_test')    c('JSPH_train','TCGA_LGG','TCGA_GBM','RBRT')
df_orig <- filter(df_orig, cohort %in% c('JSPH_test')
)

# 筛选特征
names(df_orig)
X <- colnames(df_orig)[c(1:23)] #例：这里选择了3-10号变量
df_orig <- df_orig[, which(names(df_orig)  %in% X)]
# df_orig <- na.omit(df_orig)

# 查看数据数据性质
str(df_orig)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 处理多分类,创建虚拟变量
df_orig$Grade <- factor(df_orig$Grade,levels = c("WHO2", "WHO3", "WHO4"))
levels(df_orig$Grade)
dummy_Grade <- data.frame(model.matrix(~ Grade - 1, data = df_orig))
dummy_Grade$GradeWHO4 <- NULL    # 要把哪一列作为基准列，如果这里不删除的话，默认是以因子的最后一个作为基准列
dummy_Grade$GradeWHO2 <- factor(dummy_Grade$GradeWHO2,levels = c(0, 1))
dummy_Grade$GradeWHO3 <- factor(dummy_Grade$GradeWHO3,levels = c(0, 1))
colnames(dummy_Grade)


df_orig$co_del_1p_19q <- factor(df_orig$co_del_1p_19q,levels = c("cod", "non-cod",'not_available'))
dummy_1p_19q <- data.frame(model.matrix(~ co_del_1p_19q - 1, data = df_orig))
dummy_1p_19q$co_del_1p_19qnon.cod <- NULL    # 要把哪一列作为基准列，如果这里不删除的话，默认是以因子的最后一个作为基准列
dummy_1p_19q$co_del_1p_19qcod <- factor(dummy_1p_19q$co_del_1p_19qcod,levels = c(0, 1))
dummy_1p_19q$co_del_1p_19qnot_available <- factor(dummy_1p_19q$co_del_1p_19qnot_available,levels = c(0, 1))
colnames(dummy_1p_19q)


df_orig$Histology <- factor(df_orig$Histology,levels = c("Astrocytoma", "Oligodendroglioma", "Oligoastrocytoma",'Glioblastoma'))
dummy_Histology <- data.frame(model.matrix(~ Histology - 1, data = df_orig))
dummy_Histology$HistologyGlioblastoma <- NULL    # 要把哪一列作为基准列，如果这里不删除的话，默认是以因子的最后一个作为基准列
dummy_Histology$HistologyAstrocytoma <- factor(dummy_Histology$HistologyAstrocytoma,levels = c(0, 1))
dummy_Histology$HistologyOligodendroglioma <- factor(dummy_Histology$HistologyOligodendroglioma,levels = c(0, 1))
dummy_Histology$HistologyOligoastrocytoma <- factor(dummy_Histology$HistologyOligoastrocytoma,levels = c(0, 1))
colnames(dummy_Histology)


df_orig$gender <- factor(df_orig$gender,levels = c("male", "female"))

df_orig$IDH <- factor(df_orig$IDH,levels = c("Mutant","Wild-type"))



# 其他要转换的因子变量
factor_car = c('LR','centerline','l_frontal','l_temproal','l_parietal','l_occipital','l_Insular','l_other')
for (fac in factor_car){
  df_orig[,fac] <- factor(df_orig[,fac],levels = c(0, 1))
}

cox_data <- cbind(df_orig, dummy_Grade
                  ,dummy_1p_19q
                  ,dummy_Histology
                  )
str(cox_data)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# ████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 单因素回归
names(cox_data)
X=c("RiskScore","gender",'age','IDH'
    # ,'co_del_1p_19q'
    # ,'level'
    # ,'Grade'
    ,'LR','centerline','l_frontal','l_temproal','l_parietal','l_occipital','l_Insular','l_other'
    ,'GradeWHO2','GradeWHO3'
    ,'co_del_1p_19qcod','co_del_1p_19qnot_available'
    ,'HistologyAstrocytoma','HistologyOligodendroglioma','HistologyOligoastrocytoma'
    )

# coxph(Surv(time, status) ~ Grade, data = df_orig)

univ_formulas <- sapply(X,function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = cox_data)})

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
univ_resTable
write.csv(univ_resTable,'→→单因素回归结果.csv')

# 删除空行
univ_resTable <- na.omit(univ_resTable,quote = F)


# 筛选掉 p>0.05 的特征
univ_resTable[,'p.value'] = as.numeric(univ_resTable[,'p.value']) # 修改 p.value 列的数据类型
univ_resTable <- univ_resTable[-which(univ_resTable$p.value > 0.05),]
univ_resTable

# ████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 多因素Cox回归
all_var <- row.names(univ_resTable)
all_var

# # 要删除的元素
# # RiskScore  RiskScore_log  RiskScore_17  RiskScore_17_log
# to_remove <- c('IDH'
#                # ,'co_del_1p_19q'
#                # ,'l_temproal','l_parietal','l_occipital','l_Insular','l_other'
#                )
# # 删除元素
# all_var <- all_var[!all_var %in% to_remove]
# all_var


# all_var=c("RiskScore",'age','IDH'
#     ,'co_del_1p_19q'
#     # ,'level'
#     ,'Grade'
#     # ,'LR','centerline','l_Insular','l_other'
#     ,'l_frontal','l_temproal','l_parietal','l_occipital'
#     # ,'GradeWHO2','GradeWHO3'
#     # ,'co_del_1p_19qcod','co_del_1p_19qnot_available'
#     # ,'HistologyAstrocytoma','HistologyOligodendroglioma','HistologyOligoastrocytoma'
# )

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 创建df用于多因素cox回归
muticox_data <- cox_data[,c('time','status',all_var)]
dim(muticox_data)


mul_cox <- coxph(Surv(time, status) ~ ., data = muticox_data)
# mul_cox

#一-2 multi1：提取：变量+HR+95%CI+95%CI
mul_cox1 <- summary(mul_cox)
colnames(mul_cox1$conf.int)
# multi1<-as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
multi1<-as.data.frame(mul_cox1$conf.int[, c(1, 3, 4)])
#一-3、multi2：提取：HR(95%CI)和P
multi2<-ShowRegTable(mul_cox, 
                     exp=TRUE,
                     digits=2,
                     pDigits =3,
                     printToggle = TRUE,
                     quote=FALSE, 
                     ciFun=confint)
#一-4.将两次提取结果合并成表；取名result
result <-cbind(multi1,multi2);result
#一-5.行名转为表格第一列，并给予命名"Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
View(result)

write.csv(result,'→→多因素回归结果.csv')
# p值在这儿
mul_cox

