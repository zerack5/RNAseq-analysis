
install.packages("DESeq2")

BiocManager::install("DESeq2")


path0 <- "D:/数据分析中间文件/il17a"
setwd(path0)


file1 <- "1phase_il17acount.csv"
file2 <- "2phase_il17acount.csv"
file3 <- "3phase_il17acount.csv"


library(DESeq2)
library(dplyr)
#导入counts数据矩阵，以行为基因，列为样本


dataa1 <- read.csv(file1, sep = ',',header = T,row.names = 1)

line1 <- data.frame(colnames(dataa1))
colnames(line1) <- "样本"

dataa2 <- read.csv(file2, sep = ',',header = T,row.names = 1)

line2 <- data.frame(colnames(dataa2))
colnames(line2) <- "样本"

dataa3 <- read.csv(file3, sep = ',',header = T,row.names = 1)

line3 <- data.frame(colnames(dataa3))
colnames(line3) <- "样本"

linef <- rbind(line1,line2,line3)
linef$样本

write.csv(linef,"注释信息a.csv",quote = F,row.names = T,col.names = T)

#################################################################################
data1 <- read.csv('1phase_il17a.csv', sep = ',',header = T,row.names = 1)
data2 <- read.csv('2phase_il17a.csv', sep = ',',header = T,row.names = 1)
data3 <- read.csv('3phase_il17a.csv', sep = ',',header = T,row.names = 1)

cbadata <- cbind(data1,data2)

count <- cbadata
count <- round(count)
## 过滤在所有重复样本中小于1的基因，表达量太低也没研究意义
#count <- count[rowMeans(count)>1,]
##载入样本信息
fzdata2 <- read.csv('il17a注释信息.csv', sep = ',',header = T,row.names = 1)

#fzdata2 <- read.table("C:/Users/a/Desktop/表达矩阵/注释信息.txt",header = T,row.names = 1)

#一定要变为因子数据，否者用DESeq2包分析时候会出错
fzdata2[,1] <- as.factor(fzdata2$type)

fazdata1 <- subset(fzdata2,type==1)
fazdata2 <- subset(fzdata2,type==2)
fazdata12 <- rbind(fazdata1,fazdata2)


all(rownames(fazdata12) %in% colnames(count))
all(rownames(fazdata12) == colnames(count))

dds <-  DESeqDataSetFromMatrix(countData = count,colData = fazdata12,design = ~type)
dim(dds)

#过滤
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds) 

## 差异比较
dep <- DESeq(dds)
res <- results(dep)
diff = res
diff <- na.omit(diff)  ## 去除缺失值NA
dim(diff)


write.csv(diff,"all_diff12.csv")

#Padj是P值矫正之后的数值，一般选取小于等于0.05（显著差异）的基因；同时log2FC是基因表达量的差异倍数。例如log2FC为1，证明这个基因在两种不同处理中的表达量相差了一倍，通常以大于1或小于-1为标准，大于1的为上调表达，少于-1的为下调表达。
foldChange = 1
padj = 0.05
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
dim(diffsig)

write.csv(diff,"all_diffsig12.csv")



