# 读写nifty数据有两种包来做，一个是 oro.nifti结合RNifti包，另一个是neuroim2包
# 后者安装更麻烦，但好在图像方向是对的，而且数据无误
# 前者问题在于方向和保存的问题

# install.packages('oro.nifti')
# install.packages('WhiteStripe')

# install.packages('installr')
# devtools::install_github('talgalili/installr')
# library(installr)
# install.Rtools()
# https://bbuchsbaum.github.io/neuroim2/articles/index.html
# install.packages("devtools")
# install.packages("rlang")
# update.packages(ask = FALSE)
devtools::install_github("bbuchsbaum/neuroim2")

rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(RNifti)
library(oro.nifti)
library(WhiteStripe)
library(neuroim2)
# ??RNifti
# ??WhiteStripe
# ??neuroim2



# ████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 单个WhiteStripe
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ████████████████████████████████████████████████████████████████████████
# # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 使用 oro.nifti RNifti 包读写
# 读取文件和展示图
img <- readNIfTI('t2_N4.nii.gz')
orthographic(img) # 展示图
# 
# # ——————————————————————————————————————————————
# 查看原始直方图
vals = img[img > 0]
hist(vals, breaks = 2000)
# 
# # ——————————————————————————————————————————————
# # 处理
ws = whitestripe(img = img, type = "MD", stripped = TRUE)
# names(ws)
norm = whitestripe_norm(img = img, indices = ws$whitestripe.ind)
# ?whitestripe
# # ——————————————————————————————————————————————
# 查看原始直方图 - 加了白质划线
hist(vals, breaks = 2000)
abline(v = ws$mu.whitestripe, col = "blue")
abline(v = ws$whitestripe, col = "red")

# 查看mask
mask = ws$mask.img
mask[mask == 0] = NA
orthographic(x = img, y = mask, col.y = "red")

# 查看处理后的直方图
norm_vals = norm[img > 0]
hist(norm_vals, breaks = 2000)
# 
# # 保存
# writeNifti(norm, 'demo-WSed.nii.gz', template = NULL, datatype = "double", version = 1)
# # ?writeNifti


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 使用neuroim2包
vol <- read_vol('t2_N4.nii.gz')
orthographic(vol) # 展示图，这里跟readNIfTI读出来的又是反的了

# 查看原始直方图
vals = vol[vol > 0]
hist(vals, breaks = 2000)

# ?whitestripe
# 处理
ws_neuroim2 = whitestripe(img = vol, type = "MD", stripped = TRUE)
# names(ws)
norm_neuroim2 = whitestripe_norm(img = vol, indices = ws_neuroim2$whitestripe.ind)

# 查看处理后的直方图
norm_vals = norm_neuroim2[vol > 0]
hist(norm_vals, breaks = 2000)

write_vol(norm_neuroim2, "t2_N4-WSed.nii.gz")


# ████████████████████████████████████████████████████████████████████████
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 批量WhiteStripe，使用neuroim2包
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ████████████████████████████████████████████████████████████████████████

df <- read.table('→目标文件列表-tcgaLgg-N4ed2.txt',quote = "",sep = "")   # 待处理特征

for (line in 1:nrow(df)){
  file_name <- paste0(df[line,1])
  nick_name = strsplit(file_name, '\\\\')[[1]][4]
  save_name = paste0(substr(file_name,0,nchar(file_name)-7),'-WSed.nii.gz')
  cat(file_name,'\n')
  cat(nick_name,'\n')
  
  # img = oro.nifti::readNIfTI(file_name, reorient = FALSE);
  
  # # 读取文件和展示图
  img <- read_vol(file_name)
  # 处理
  ws = whitestripe(img = img, type = "MD", stripped = TRUE)
  norm = whitestripe_norm(img = img, indices = ws$whitestripe.ind)
  # 保存
  write_vol(norm, save_name)
  cat('whiteStripe succeed!\n')

  rm(file_name,nick_name,save_name,img,ws,norm)
}




