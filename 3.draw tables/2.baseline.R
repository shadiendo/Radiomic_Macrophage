# install.packages('compareGroups')

rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(table1) 
library(compareGroups)


# █████████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 基础三线表绘制
# http://t.csdn.cn/cjXjb
# http://t.csdn.cn/prY4g
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
df_orig <- read.csv('GeneralView.csv', header=T,row.names = 1, encoding="UTF-8")   # 影像组学数据
names(df_orig)
# 先因子化
factor_demo = c('gender','IDH','status','co_del_1p_19q','Histology')
df_orig[,factor_demo]<-lapply(df_orig[,factor_demo],as.factor)


# 理顺y轴方向，逆向则应该对levels和labels： rev(content_long$GO_ID)
df_orig$cohort <- factor(df_orig$cohort,
                         levels = df_orig$cohort,
                         labels = df_orig$cohort)

# # 设置年龄
# age_class <- c()
# for (i in df_orig$age){
#   if (i == 'Not_Available'){
#     age_class <- c(age_class,'Not_Available')
#   }else if(i <45){
#       age_class <- c(age_class,"<45")
#   }else if(i >=45){
#         age_class <- c(age_class,'≥45')}
# }
# df_orig$age_class <- age_class




# 设置因子的表示
df_orig$cohort <- factor(df_orig$cohort,
    levels = c('JSPH_train','TCGA','RBRT','JSPH_test'),
    labels = c("JSPH(Discovery)","TCGA(Discovery)","RBRT(Discovery)","JSPH(Validation)"))
df_orig$status <- factor(df_orig$status,
    levels = c('0','1'),
    labels = c("Alive","Dead"))
df_orig$gender <- factor(df_orig$gender,
    levels = c('male','female'),
    labels = c("Male","Female"))
df_orig$Histology <- factor(df_orig$Histology,
    levels = c('Astrocytoma','Oligodendroglioma','Oligoastrocytoma','Glioblastoma'),
    labels = c('Astrocytoma','Oligodendroglioma','Oligoastrocytoma','Glioblastoma'))
df_orig$IDH <- factor(df_orig$IDH,
  levels = c('Mutant','Wild-type'),
  labels = c("Mutation","Wild-type"))
df_orig$co_del_1p_19q <- factor(df_orig$co_del_1p_19q,
    levels = c('cod','non-cod','--'),
    labels = c("Co-deletion","Non Co-deletion","NA"))

summary(df_orig)


# 设置大标签，这里的标签就是行名。否则三线表行名会显示为status，age等等
label(df_orig$age_abs) <- "Age"
label(df_orig$gender) <- "Sex"
# label(df_orig$coDelete_1p19q) <- "1p19q status"
label(df_orig$IDH) <- "IDH1 status"
label(df_orig$status) <- "Status"
label(df_orig$time) <- "Time"
label(df_orig$co_del_1p_19q) <- "1p19q status"

#设置单位，这样可以在行名那儿显示单位，如图中的（years，Kg）
units(df_orig$age_abs) <- 'year'



# 多层次绘图
table1(~ age_abs + gender + status + time+
         Histology + Grade +
         IDH + co_del_1p_19q+
         LR+centerline + l_frontal+l_temproal+l_parietal+l_occipital+l_Insular+l_other
       | cohort, data=df_orig,
       # overall = '【TOTAL】',   # 是否关闭汇总
       overall = F,   # 是否关闭汇总
       # topclass="Rtable1-zebra",
       topclass="Rtable1-zebra Rtable1-shade"
       )



# # 使用descrTable函数
# formula = cohort ~ age_abs + gender + status + time+ 
#   Histology + Grade + 
#   IDH + co_del_1p_19q+
#   LR+centerline + l_frontal+l_temproal+l_parietal+l_occipital+l_Insular+l_other
# 
# descrTable(formula = formula,
#            # show.all = T,
#            show.p.ratio = T,
#            data = df_orig,
# )



# 使用descrTable函数
formula = class ~ age_abs + gender + status + time+
  Histology + Grade +
  IDH + co_del_1p_19q+
  LR+centerline + l_frontal+l_temproal+l_parietal+l_occipital+l_Insular+l_other

pppp <- descrTable(formula = formula,
           # show.all = T,
           show.p.ratio = T,
           data = df_orig,
)
export2csv(pppp, file='baseline_p.csv')








