rm(list=ls())
# install.packages('pheatmap')
# install.packages('ggdendro')
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(pheatmap)
library(dplyr)
library(tidyverse)
library(ggdendro)


color_pool = colorRampPalette(rev(c(  "#ed4465",  # -1
                                      '#ed4f6b',  # -0.2
                                      '#ef5e72',  # -0.4
                                      '#f06e7a',  # -0.6
                                      '#f38084',  # -0.8
                                      "#f4908c",  # 0  中间
                                      '#f7a396',  # -0.2
                                      '#f9b59f',  # -0.4
                                      '#fac6a8',  # -0.5
                                      '#fcd6b0',  # -0.6
                                      "#fee7b9"  # 1
)))(50)



# ████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 读取数据
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# open_file = "1000iCoxRet.csv"
open_file = "1000iCoxLasso.csv"
# open_file = "1000iRSFVH.csv"

file_theme = substr(open_file,1,nchar(open_file)-4)
content_orig <- read.table(file=open_file, header = TRUE, sep = ",",row.names=1)


# 筛选特别低的
feature_preselect_ord=c()
threshold_least = 20
for (row in 1:nrow(content_orig)){
  flag=FALSE        # 要不要把当前特征选择
  for (col in 1:ncol(content_orig)){
    cell = content_orig[row,col]
    if (cell>threshold_least){
      flag=TRUE
    }
  }
  if (flag==TRUE){feature_preselect_ord = c(feature_preselect_ord,row)}
}

content <- content_orig[feature_preselect_ord,1:ncol(content_orig)] %>%
  scale()%>%    #缩放数据
  t()

content_orig <- content_orig %>%
  scale()%>%
  t()

rm(cell,col,flag,row,threshold_least,feature_preselect_ord)

# ████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 运行热图筛选
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
mark_the_heatpamt <- function(c_method, cut_num, title_demo, group_mark){
  # c_method='ward.D2'
  # cut_num=5
  # title_demo="CoxML"
  # group_mark = c(1,3,5)

  # complete(默认)、mcquitty、ward.D、ward.D2、single、median、centroid、average
  p<-pheatmap(content,
              main = title_demo,
              color = color_pool,
              show_rownames = TRUE,show_colnames = FALSE,   # 行列名称
              cluster_rows = FALSE,cluster_cols = TRUE,   # 行列聚类
              clustering_method = c_method,
              
              # legend = FALSE,
              treeheight_col=30,  # 行集群树的高度
              
              cutree_cols = cut_num, # 切割列
              border_color = NA,  # 不显示网格线
  )
  p
  
  # ————————————————————————————————————————
  # 聚类切割情况
  row_cluster <- as.data.frame(cutree(p$tree_col, k=cut_num))  # 获取各特征的聚类分类，顺序是原顺序
  row_cluster = rownames_to_column(row_cluster,var='name')  # 使用dplyr和tidyverse包的函数提取行名
  # 按聚类顺序排个序
  row_cluster <- as.data.frame(row_cluster[p$tree_col$order,],row.names = rownames(row_cluster)) 
  colnames(row_cluster)[2]="Cluster"   # 重命名分组的标题
  row_cluster
  
  # ——————————————————————————————————————————————
  # 查看各组情况，及为再次绘图做准备
  gatherShow = as.data.frame(table(row_cluster$Cluster))
  gatherShow = gatherShow[unique(row_cluster$Cluster),]
  gatherShow
  
  
  # 对分组汇总情况添加信息
  for (i in gatherShow$Var1){
    gatherShow[i,3] = paste('n=',as.character(gatherShow[i,2]),sep = "", collapse = '')
  }
  colnames(gatherShow)[3]="Cluster_name" 
  
  # 设置颜色信息，标记重点颜色
  keycolor = 'red4'
  gatherShow$col_color = rep('gray90',nrow(gatherShow))
  # rownames(gatherShow) <- 1:nrow(gatherShow)
  gatherShow$col_color[group_mark] <- keycolor  # 根据要求来标记
  gatherShow
  
  # 设置特征的分组信息col_df
  col_df <- as.data.frame(cutree(p$tree_col, k=cut_num))
  colnames(col_df)[1]="Cluster" 
  # 将原本col_df中表示分组的类似1，转换从n=117这种形式
  for (i in gatherShow$Var1){col_df[col_df == i] = gatherShow[i,3]}; rm(i);
  # View(col_df)
  
  # 设置分组的颜色
  Cluster <- gatherShow$col_color
  names(Cluster) <- gatherShow$Cluster_name
  ann_colors = list(Cluster)
  names(ann_colors) = 'Cluster'
  ann_colors
  
  # 绘图
  p<-pheatmap(content, 
              main = title_demo,
              color = color_pool,
              show_rownames = TRUE,show_colnames = FALSE,   # 行列名称
              cluster_rows = FALSE,cluster_cols = TRUE,   # 行列聚类
              # legend = FALSE,
              treeheight_col=30,  # 行集群树的高度
              
              clustering_method = c_method,
              
              cutree_cols = cut_num,  # 切割列
              
              annotation_col = col_df,
              
              annotation_colors = ann_colors,
              
              border_color = NA
  )
  p
  
  # 提取出哪个类的特征，并保存
  group_id_selected = gatherShow[group_mark,]$Var1
  features_by_hotmap = filter(row_cluster,row_cluster$Cluster %in% group_id_selected)[,1]

  # 保存那行特征名
  features_by_hotmap_save = paste('→',title_demo,"_Selected_num",length(features_by_hotmap),'.txt', sep = "")
  write.table(as.data.frame(features_by_hotmap), features_by_hotmap_save,row.names =F, col.names=F,quote =F)
  cat('选出特征',length(features_by_hotmap),'个','保存至【',features_by_hotmap_save,'】')
  return(p)
}

# pdf保存函数
save_pheatmap_pdf <- function(filename,plot, width, height) {
  stopifnot(!missing(plot))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(plot$gtable)
  dev.off()
}

# complete(默认)、mcquitty、ward.D、ward.D2、single、median、centroid、average
plotHeatmap_CoxML = mark_the_heatpamt('ward.D2',4,file_theme,c(1))


save_pheatmap_pdf(paste0('→→',file_theme,'_Heatmap.pdf'),plotHeatmap_CoxML,7,4)


# ████████████████████████████████████████████████████████████████████████
# rm(list=ls())
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 特征降维图
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# install.packages('umap')
# install.packages('Rtsne')
# install.packages('cowplot')
library(ggplot2)
library(umap)
library(Rtsne)
library(cowplot)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1.确定哪些特征待处理
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
rf_selected <- read.csv('→1000iCoxLasso_Selected_num35.txt', sep='\t',header=F)$V1

DimensionalityReduction_plot <- function(whitchEntirety,rf_selected_step1,title_demo){
  # whitchEntirety = content_orig
  # rf_selected_step1 <- rf_selected
  # title_demo = 'CoxML'
  
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # tsne
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  df_orig <- as.data.frame(t(whitchEntirety))
  df_orig = unique(df_orig) # 去除重复数据

  # 分组
  df_orig$Group <- 'ABANDON'
  for (i in rf_selected_step1){df_orig[i,'Group'] <- 'SELECTED'}

  # 删除有空值的行
  df_orig <- na.omit(df_orig)

  # 开始tSNE
  estimate_df = df_orig
  estimate_df = estimate_df[,1:length(estimate_df)-1]

  # set.seed(111) # 设置随机数种子，此处与PCA不同，t-SNE的结果具有随机性，所以为了每次结果都可以重复就要设置随机种子
  tsne_ret <- Rtsne(estimate_df,
               dims = 2, # 维度2，即x和y轴
               max_iter = 1000,# max_iter 为迭代次数，默认为1000，因为此数值过小会导致结果未完全收敛，过大会浪费计算量，导致等待时间变长。
               theta = 0.5, # 计算速度与精确度之间的权衡，范围在0~1之间，越接近0越精确，默认0.5。不同的值最终出来的结果也不同，大家可以自己改变一下看看最后结果怎么样，选取自己最喜欢的结果。
               perplexity = 20, # 满足 3*perplexity < nrow(data) - 1 ,随着 perplexity 值的增加，形状越来越清晰。
               verbose = T)# 运行大量数据的时候可以看计算到了什么程度


  # 将画图需要的数据从上一步整理出来，并调整好内容
  site <- data.frame(tsne_ret$Y) #将数据中的x y轴中的数据抽取出来
  colnames(site) <- c("tSNE1","tSNE2")#将site的列名更改

  # 分组
  site$Group <- df_orig$Group


  # 作图
  p <- ggplot(site, aes(tSNE1, tSNE2, fill = Group, shape = Group)) +
    ggtitle(title_demo) +
    geom_point(size = 2,alpha=1)+
    #添加置信椭圆，注意不是聚类
    stat_ellipse(aes(color = Group), geom = 'polygon', level = 0.95, alpha = 0.05, linetype = 2, show.legend = TRUE) +

    # 填充的颜色
    scale_fill_manual(values = c("#00a1d5","#b24745","violet"))+
    # 边的颜色
    scale_color_manual(values = c("#00a1d5","#b24745","violet"))+
    # 不同数字代表不同的点样式
    scale_shape_manual(values = c(21,22,24))+


    # 图例位置"left" "right" "bottom" "top" 或者自己定义 c(0.9,0.7)
    # theme(legend.position = c(0.9,0.1))+
    theme(legend.position = 'bottom')+

    # 网格是否要去除
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))+
    # 背景颜色
    theme_bw()
  p


  return(p)
  
  

  # # # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # # # umap
  # # # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # df_orig <- as.data.frame(t(whitchEntirety))
  # 
  # # 分组
  # df_orig$Group <- 'ABANDON'
  # for (i in rf_selected_step1){df_orig[i,'Group'] <- 'SELECTED'}
  # 
  # estimate_df = df_orig
  # estimate_df = estimate_df[,1:length(estimate_df)-1]
  # 
  # # 正式运算
  # umap <-  umap(estimate_df)
  # # # ?umap()
  # #
  # # 将画图需要的数据从上一步整理出来，并调整好内容
  # site <- data.frame(umap$layout) #将数据中的x y轴中的数据抽取出来
  # colnames(site) <- c("UMAP1","UMAP2")#将site的列名更改
  # 
  # # 分组
  # site$Group <- df_orig$Group
  # 
  # # 使用ggplot2作图
  # p <- ggplot(site, aes(UMAP1, UMAP2, fill = Group, shape = Group)) +
  #   ggtitle(title_demo) +
  #   geom_point(size = 2,alpha=1)+#添加点，size为点的大小
  # 
  #   #添加置信椭圆，注意不是聚类
  #   stat_ellipse(aes(color = Group), geom = 'polygon',
  #                level = 0.95, alpha = 0.15, linetype = 2,
  #                show.legend = TRUE) +
  # 
  # 
  #   # 填充的颜色
  #   scale_fill_manual(values = c("#00a1d5","#b24745","violet"))+
  #   # 边的颜色
  #   scale_color_manual(values = c("#00a1d5","#b24745","violet"))+
  #   # 不同数字代表不同的点样式
  #   scale_shape_manual(values = c(21,22,24))+
  # 
  # 
  # 
  #   # 网格是否要去除
  #   theme(panel.border = element_blank(),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         axis.line = element_line(colour = "black"))+
  #   # 背景颜色
  #   theme_bw()
  # 
  #   # # 背景颜色
  #   # theme_bw()+
  #   # # 图例位置"left" "right" "bottom" "top" 或者自己定义 c(0.9,0.7)
  #   # theme(legend.position = c(0.9,0.1))+
  #   # theme(panel.border = element_blank(),
  #   #       panel.grid.major = element_blank(),
  #   #       panel.grid.minor = element_blank(),
  #   #       axis.line = element_line(colour = "black"))
  # p
  # return(p)

}

# 这里选择的整体应该为 content_orig
plot_tSNE_CoxML = DimensionalityReduction_plot(content_orig,rf_selected,file_theme)
plot_tSNE_CoxML

ggsave(paste0('→→→',file_theme,'_tSNE.pdf'), plot_tSNE_CoxML, width =9, height = 4, dpi = 300)



# █████████████████████████████████████████████████████████████████████████████

# 核密度图

# █████████████████████████████████████████████████████████████████████████████
rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1.数据准备
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
file_name <- '1000iRSFVH.csv'
df_orig <- read.csv(file_name, header=T, row.names=1)   # 打开文件  , row.names=1    , fileEncoding = 'GB2312'
colnames(df_orig)
# fix(df_orig)

# 筛选特别低的
feature_preselect_ord=c()
threshold_least = 50
for (row in 1:nrow(df_orig)){
  flag=FALSE        # 要不要把当前特征选择
  for (col in 1:ncol(df_orig)){
    cell = df_orig[row,col]
    if (cell>threshold_least){
      flag=TRUE
    }
  }
  if (flag==TRUE){feature_preselect_ord = c(feature_preselect_ord,row)}
}
df_orig <- df_orig[feature_preselect_ord,1:ncol(df_orig)]


df_glioma <- df_orig%>%dplyr::select(c('parm37','parm46','parm55','parm64','parm73'))
df_lgg <- df_orig%>%dplyr::select(c('parm37_lgg','parm46_lgg','parm55_lgg','parm64_lgg','parm73_lgg'))
df_gbm <- df_orig%>%dplyr::select(c('parm37_gbm','parm46_gbm','parm55_gbm','parm64_gbm','parm73_gbm'))
df_con <- df_orig%>%dplyr::select(c('parm37_conjugate','parm46_conjugate','parm55_conjugate','parm64_conjugate','parm73_conjugate'))



# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2.绘图函数
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
showPlot<- function(df_demo,titleName_demo,haveLegend){
  # df_demo = df_glioma
  # titleName_demo='lgg&gbm'
  # haveLegend = 'openLegend'
  
  
  # 转化为ggplot2格式的数据类型，长数据
  df_demo_long = df_demo %>%
    rownames_to_column("RF_names") %>%
    pivot_longer(-1,names_to = "class",values_to = "Value")
  
  # # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # # 绘制核密度图 ggplot 单线
  # # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # ggplot(df_demo) +
  #   geom_density(aes(perm37,y=..density..),size=1,color="black")+
  #   # theme_minimal()+
  #   theme(panel.grid.major=element_blank())
  
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # 绘制核密度图 ggplot 多线
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  gg <- ggplot(df_demo_long) +
    ggtitle(titleName_demo)+
    xlab("Number Of Repeats") +
    ylab("Kernel Density") +
    geom_density(alpha=.15,
                 aes(x=Value,
                     # color = class,  # 线条颜色
                     fill = class    # 填充的颜色
                 )
    )+
    # coord_cartesian(ylim = c(0,0.05),expand = TRUE)+     # 缩放y轴局部
    
    theme(plot.title = element_text(color="black", size=18, face="bold",hjust=0.5),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.y = element_text(size = 12,colour="black"),
          axis.text.x = element_text(size = 12,colour="black"))
  
  if (haveLegend=='closeLegend'){
    gg = gg+scale_fill_discrete(guide=FALSE)
  }else if (haveLegend=='openLegend'){
    gg = gg+ scale_fill_discrete(name="class",
                                 labels=c("perm37", "perm46", "perm55", "perm64", "perm73"))
  }
  return(gg)
  
}


p1 = showPlot(df_glioma,'glioma','closeLegend')
p2 = showPlot(df_lgg,'lgg','closeLegend')
p3 = showPlot(df_gbm,'gbm','closeLegend')
p4 = showPlot(df_con,'conjugate','openLegend')


ret = p1|p2|p3|p4
ret


ggsave(paste0('RSFVH','_HeMiDu.pdf'), ret, width =24, height = 3, dpi = 300)












