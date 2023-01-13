library(plyr)



aimp <- "SOX9"


setwd("C:/Users/a/Desktop/表达矩阵")
b3gnt2as <- read.csv('exp2as.csv', sep = ',',header = T,row.names = 1)
b3gnt2c <- read.csv('exp2control.csv', sep = ',',header = T,row.names = 1)
nkb3gnt2as <- b3gnt2as[c(grep("_N",colnames(b3gnt2as)))]

a2nk <- nkb3gnt2as[aimp,]

nkb3gnt2c <- b3gnt2c[c(grep("_N",colnames(b3gnt2c)))]
b2nk <- nkb3gnt2c[aimp,]





setwd("C:/Users/a/Desktop/表达矩阵")
b3gnt2as <- read.csv('exp2as.csv', sep = ',',header = T,row.names = 1)
b3gnt2c <- read.csv('exp2control.csv', sep = ',',header = T,row.names = 1)
nkb3gnt2as <- b3gnt2as[c(grep("_M",colnames(b3gnt2as)))]

a2mo <- nkb3gnt2as[aimp,]

nkb3gnt2c <- b3gnt2c[c(grep("_M",colnames(b3gnt2c)))]
b2mo <- nkb3gnt2c[aimp,]





setwd("C:/Users/a/Desktop/表达矩阵")
b3gnt2as <- read.csv('exp2as.csv', sep = ',',header = T,row.names = 1)
b3gnt2c <- read.csv('exp2control.csv', sep = ',',header = T,row.names = 1)
nkb3gnt2as <- b3gnt2as[c(grep("CD4",colnames(b3gnt2as)))]

a2cd4 <- nkb3gnt2as[aimp,]

nkb3gnt2c <- b3gnt2c[c(grep("CD4",colnames(b3gnt2c)))]
b2cd4 <- nkb3gnt2c[aimp,]



setwd("C:/Users/a/Desktop/表达矩阵")
b3gnt2as <- read.csv('exp2as.csv', sep = ',',header = T,row.names = 1)
b3gnt2c <- read.csv('exp2control.csv', sep = ',',header = T,row.names = 1)
nkb3gnt2as <- b3gnt2as[c(grep("CD8",colnames(b3gnt2as)))]

a2cd8 <- nkb3gnt2as[aimp,]

nkb3gnt2c <- b3gnt2c[c(grep("CD8",colnames(b3gnt2c)))]
b2cd8 <- nkb3gnt2c[aimp,]


cd8 <- rbind.fill(a2cd8,b2cd8)
cd4 <- rbind.fill(a2cd4,b2cd4)
cd <- rbind.fill(cd8,cd4)
mo <- rbind.fill(a2mo,b2mo)
nk <- rbind.fill(a2nk,b2nk)

allsort <- rbind.fill(cd,mo,nk)
row.names(allsort) <- c("as cd8","control cd8","as cd4","control cd4","as mo","control mo","as nk","control nk")


#####画小提琴图

library(ggplot2)
library(reshape2)
library(ggsignif)

# 读取小提琴图数据文件
df <- allsort
# 把数据转换成ggplot常用的类型（长数据）
df <- t(df)

df = melt(df)

library(tidyverse)
df <- drop_na(df)


# cd4,cd8细胞绘图
dfcd <- df[c(grep(pattern="cd",df$Var2)),]


ggplot(dfcd,aes(x=Var2,y=value,fill=Var2))+
  geom_violin(alpha = 1,              # 透明度
              trim = T,               # 是否修剪尾巴，即将数据控制到真实的数据范围内
              scale = "area",         # 如果“area”(默认)，所有小提琴都有相同的面积(在修剪尾巴之前)。如果是“count”，区域与观测的数量成比例。如果是“width”，所有的小提琴都有相同的最大宽度。
              
  )+
  theme_bw()+                          # 白色主题
  theme(
    axis.text.x = element_text(angle = 90,
                               vjust = 0.5
    )       # x轴刻度改为倾斜90度，防止名称重叠
  )+
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("as cd8","control cd8"),c("as cd4","control cd4")), # 选择你想在哪2组上添加标签
    step_increase = 0.1,
    test="t.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    test.args = list("var.equal" = T),  # 等方差    
    map_signif_level=F                 # 标签样式F为数字，T为*号
  )



#nk和mo细胞绘图





dfnkmo <- df[c(grep(pattern="mo|nk",df$Var2)),]


ggplot(dfnkmo,aes(x=Var2,y=value,fill=Var2))+
  geom_violin(alpha = 1,              # 透明度
              trim = T,               # 是否修剪尾巴，即将数据控制到真实的数据范围内
              scale = "area",         # 如果“area”(默认)，所有小提琴都有相同的面积(在修剪尾巴之前)。如果是“count”，区域与观测的数量成比例。如果是“width”，所有的小提琴都有相同的最大宽度。
              
  )+
  theme_bw()+                          # 白色主题
  theme(
    axis.text.x = element_text(angle = 90,
                               vjust = 0.5
    )       # x轴刻度改为倾斜90度，防止名称重叠
  )+
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("as mo","control mo"),c("as nk","control nk")), # 选择你想在哪2组上添加标签
    step_increase = 0.1,
    test="t.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    test.args = list("var.equal" = T),  # 等方差    
    map_signif_level=F                 # 标签样式F为数字，T为*号
  )









#总的细胞绘图

ggplot(df,aes(x=Var2,y=value,fill=Var2))+
  geom_violin(alpha = 1,              # 透明度
              trim = T,               # 是否修剪尾巴，即将数据控制到真实的数据范围内
              scale = "area",         # 如果“area”(默认)，所有小提琴都有相同的面积(在修剪尾巴之前)。如果是“count”，区域与观测的数量成比例。如果是“width”，所有的小提琴都有相同的最大宽度。
              
  )+
  theme_bw()+                          # 白色主题
  theme(
    axis.text.x = element_text(angle = 90,
                               vjust = 0.5
    )       # x轴刻度改为倾斜90度，防止名称重叠
  )+
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("as cd8","control cd8"),c("as cd4","control cd4"),c("as mo","control mo"),c("as nk","control nk")), # 选择你想在哪2组上添加标签
    step_increase = 0.1,
    test="t.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    test.args = list("var.equal" = T),  # 等方差    
    map_signif_level=F                 # 标签样式F为数字，T为*号
  )










ma2cd8 <- data.frame(matrix(ncol = 2, nrow = 1))

ma2cd8 <- []
ma2cd8$counts <- mean(a2cd8[,1])

ma2cd8$label <- aimp
ma2cd8 <- as.matrix(ma2cd8)


ma2cd4 <- mean(a2cd4)


boxplot(nkb3gnt2as["ACTB",],nkb3gnt2c["B3GNT2",],main="NK细胞表达量",xlab="种类",ylab="表达量",col = c(2,3),names =c("AS","Normal"))



