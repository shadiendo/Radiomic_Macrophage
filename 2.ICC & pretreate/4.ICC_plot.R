library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd())
library(ggplot2)

rm(list=ls())

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1.数据准备
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
file_name <- '→→→ICC2k绘图用.csv'
df_orig <- read.csv(file_name, header=T, row.names=2, fileEncoding = 'GB2312')   # 打开文件

# 把序号那列改名为 id
fix(df_orig)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 绘制statistics的散点图
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# R语言自己的包
pdf("scatter_plot.pdf", width=6, height=5)
# 设置坐标 x 轴范围 2.5 到 5, y 轴范围 15 到 30.
plot(x = df_orig$id,y = df_orig$ICC,
     xlab = "p_value",
     ylab = "percent",
     # xlim = c(2.5,5),
     # ylim = c(15,30),
     main = "Weight vs Milage"
)
dev.off()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ggplot包散点图
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 最简单的散点图
ggplot(df_orig, aes(x = id, y = ICC)) +
  # 散点图函数
  geom_point()



# 映射离散型型变量
plot_ret<- ggplot(df_orig, aes(x = id, y = ICC, color = delete)) +
  guides(color=FALSE)+   # 不显示图例
  #———————————————————————————————————————————— X轴标签
  geom_point(size = 0.5,alpha = 1)+   # 散点图函数 alpha设置散点透明度
  geom_hline(yintercept=0.9, linetype = "solid", size = 0.5)+    # 增加一条刻度线
  
  #———————————————————————————————————————————— 标题、坐标轴
  ggtitle("Screened By ICC")+

  xlab("Default Pyradiomics Extractor order") +
  scale_x_discrete(expand = waiver())+  # 取消内容与x轴左右的间距(间隙)

  ylab("Default Pyradiomics Extractor order") +
  scale_y_continuous(breaks= c(seq(0, 1, by = 0.25),0.9))+     # 添加某局部标签0.9
  coord_cartesian(ylim = c(0.03,0.96),expand = TRUE)+      # 缩放y轴局部


  theme(plot.title = element_text(color="black", size=22, face="bold",hjust=0.5),
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        axis.text.y = element_text(size = 15,colour="black"))+

  # #————————————————————————————————————————————
  scale_color_manual(values = c("black","grey"))+    # 手动设置颜色

  #ggplot2去掉网格线和背景色
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
plot_ret

ggsave("ICC_plot.pdf", plot_ret, width = 5, height = 4.8,dpi = 300)#保存文件



