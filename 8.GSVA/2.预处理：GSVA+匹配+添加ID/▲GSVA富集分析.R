# options(download.file.method="libcurl")
# options(url.method="libcurl")
# 
# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# if (!require("GSVA", quietly = TRUE)) {
#   BiocManager::install("GSVA")
# }
# if (!require("GSEABase", quietly = TRUE)) {
#   BiocManager::install("GSEABase")
# }
# if (!require(clusterProfiler)) {
#   BiocManager::install("clusterProfiler")
# }
# # 安装limma    Biobase
# if (!requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
# }
# BiocManager::install("limma")
# BiocManager::install("Biobase")
# install.packages("pheatmap")
# 
# # 下载Rtool
# https://cran.r-project.org/bin/windows/Rtools/history.html



#引用包
rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(GSVA)
library(GSEABase)
library(clusterProfiler)
library(BiocManager)
library(pheatmap)
library(limma)


exp.path ="RBRT_exp_pretreated.csv" 
gmtFile="combined.symbols.gmt"                #基因集文件

exp.data <- read.csv(exp.path,header=T, row.names=1)

data <- as.matrix(exp.data)

#GSVA分析
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
gsvaResult=gsva(data, 
                geneSets,
                method="ssgsea",
                ssgsea.norm=TRUE,
                min.sz=10, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=1)

#输出GSVA结果
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)

save_name = paste0('★GSVA_',substr(exp.path,1,nchar(exp.path)-4),'.txt')
write.table(gsvaOut, file=save_name, sep="\t", quote=F, col.names=F)



