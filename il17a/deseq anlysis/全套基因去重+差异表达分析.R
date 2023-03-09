
install.packages("DESeq2")
install.packages("limma")
BiocManager::install("DESeq2")


path0 <- "D:/数据分析中间文件/il17a"
setwd(path0)


file1 <- "1phase_il17acount.csv"
file2 <- "2phase_il17acount.csv"
file3 <- "3phase_il17acount.csv"
kname2 <- "all_diff23_paired_res.csv"
kname23 <- "all_diffsig23_paired_gene_res.csv"
kname23up <- "all_diffsig23_paired_gene_res_up.csv"
kname23down <- "all_diffsig23_paired_gene_res_down.csv"

file1 <- "1phase_il17acount.csv"
file2 <- "2phase_il17acount.csv"
file3 <- "3phase_il17acount.csv"
kname1 <- "all_diff12_paired_gene_res.csv"
kname12 <- "all_diffsig12_paired_gene_res.csv"
kname12up <- "all_diffsig12_paired_gene_res_up.csv"
kname12down <- "all_diffsig12_paired_gene_res_down.csv"






uname1 <- "g1phase_il17acount.csv"
uname2 <- "g2phase_il17acount.csv"
uname3 <- "g3phase_il17acount.csv"




library(DESeq2)
library(dplyr)
#导入counts数据矩阵，以行为基因，列为样本


dataa1 <- read.csv(file1, sep = ',',header = T,row.names = 1)

dataa1$ensembl_gene_id <-rownames(dataa1)



####id change

library(biomaRt)
listMarts()
mart<-useMart("ensembl")
dataset_list=as.data.frame(listDatasets(mart))
dataset_list
mart_oas=useMart("ensembl", "hsapiens_gene_ensembl")#human_gene

filter_list=as.data.frame(listFilters(mart_oas))
tail(filter_list)
attri_list=as.data.frame(listAttributes(mart_oas))
attri_list

list_gene=getBM(attributes = c("ensembl_gene_id","entrezgene_id","external_gene_name","description","chromosome_name"),
                filters = "ensembl_gene_id",
                values = dataa1$ensembl_gene_id,#ensemble_id
                mart = mart_oas)

t3g <- list_gene[,c("ensembl_gene_id","external_gene_name")]

datap <- merge(dataa1,t3g)



expr<- datap[,-1]

expr_mean=aggregate(.~external_gene_name,mean,data=expr)

expr_mean <- expr_mean[-1,]

rownames(expr_mean)<-expr_mean$external_gene_name
expr_mean <- expr_mean[,-1]

write.csv(expr_mean,file = uname1)
#########################################
#####file2文件转换

dataa1 <- read.csv(file2, sep = ',',header = T,row.names = 1)

dataa1$ensembl_gene_id <-rownames(dataa1)



####id change

library(biomaRt)
listMarts()
mart<-useMart("ensembl")
dataset_list=as.data.frame(listDatasets(mart))
dataset_list
mart_oas=useMart("ensembl", "hsapiens_gene_ensembl")#human_gene

filter_list=as.data.frame(listFilters(mart_oas))
tail(filter_list)
attri_list=as.data.frame(listAttributes(mart_oas))
attri_list

list_gene=getBM(attributes = c("ensembl_gene_id","entrezgene_id","external_gene_name","description","chromosome_name"),
                filters = "ensembl_gene_id",
                values = dataa1$ensembl_gene_id,#ensemble_id
                mart = mart_oas)

t3g <- list_gene[,c("ensembl_gene_id","external_gene_name")]

datap <- merge(dataa1,t3g)



expr<- datap[,-1]

expr_mean=aggregate(.~external_gene_name,mean,data=expr)

expr_mean <- expr_mean[-1,]

rownames(expr_mean)<-expr_mean$external_gene_name
expr_mean <- expr_mean[,-1]

write.csv(expr_mean,file = uname2)

#########################################
#####file3文件转换

dataa1 <- read.csv(file3, sep = ',',header = T,row.names = 1)

dataa1$ensembl_gene_id <-rownames(dataa1)



####id change

library(biomaRt)
listMarts()
mart<-useMart("ensembl")
dataset_list=as.data.frame(listDatasets(mart))
dataset_list
mart_oas=useMart("ensembl", "hsapiens_gene_ensembl")#human_gene

filter_list=as.data.frame(listFilters(mart_oas))
tail(filter_list)
attri_list=as.data.frame(listAttributes(mart_oas))
attri_list

list_gene=getBM(attributes = c("ensembl_gene_id","entrezgene_id","external_gene_name","description","chromosome_name"),
                filters = "ensembl_gene_id",
                values = dataa1$ensembl_gene_id,#ensemble_id
                mart = mart_oas)

t3g <- list_gene[,c("ensembl_gene_id","external_gene_name")]

datap <- merge(dataa1,t3g)



expr<- datap[,-1]

expr_mean=aggregate(.~external_gene_name,mean,data=expr)

expr_mean <- expr_mean[-1,]

rownames(expr_mean)<-expr_mean$external_gene_name
expr_mean <- expr_mean[,-1]

write.csv(expr_mean,file = uname3)






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

#################################################################################23
data1 <- read.csv('g1phase_il17acount.csv', sep = ',',header = T,row.names = 1)
data2 <- read.csv('g2phase_il17acount.csv', sep = ',',header = T,row.names = 1)
data3 <- read.csv('g3phase_il17acount.csv', sep = ',',header = T,row.names = 1)

cbadata <- cbind(data2,data3)

count <- cbadata
count <- round(count)
## 过滤在所有重复样本中小于1的基因，表达量太低也没研究意义
#count <- count[rowMeans(count)>1,]
##载入样本信息
fzdata2 <- read.csv('il17a注释信息.csv', sep = ',',header = T,row.names = 1)

#fzdata2 <- read.table("C:/Users/a/Desktop/表达矩阵/注释信息.txt",header = T,row.names = 1)

#一定要变为因子数据，否者用DESeq2包分析时候会出错
#fzdata2[,1] <- as.factor(fzdata2$type)

fazdata1 <- subset(fzdata2,type==2)
fazdata2 <- subset(fzdata2,type==3)
fazdata12 <- rbind(fazdata1,fazdata2)

as.data.frame(fazdata12)

all(rownames(fazdata12) %in% colnames(count))
all(rownames(fazdata12) == colnames(count))

fazdata12$type <- as.character(fazdata12$type)
fazdata12$patient_id <- as.character(fazdata12$patient_id)
dds <-  DESeqDataSetFromMatrix(countData = count,colData = fazdata12,design = ~type+patient_id)
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


write.csv(diff,kname2)

#Padj是P值矫正之后的数值，一般选取小于等于0.05（显著差异）的基因；同时log2FC是基因表达量的差异倍数。例如log2FC为1，证明这个基因在两种不同处理中的表达量相差了一倍，通常以大于1或小于-1为标准，大于1的为上调表达，少于-1的为下调表达。
foldChange = 1
padj = 0.05
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
dim(diffsig)

write.csv(diffsig,kname23)

up_DEG <- subset(diffsig, padj < 0.05 & log2FoldChange > 1)
down_DEG <- subset(diffsig, padj < 0.05 & log2FoldChange < -1)

write.csv(up_DEG,kname23up)
write.csv(down_DEG,kname23down)


#################################################################################12
data1 <- read.csv('g1phase_il17acount.csv', sep = ',',header = T,row.names = 1)
data2 <- read.csv('g2phase_il17acount.csv', sep = ',',header = T,row.names = 1)
data3 <- read.csv('g3phase_il17acount.csv', sep = ',',header = T,row.names = 1)

cbadata <- cbind(data1,data2)

count <- cbadata
count <- round(count)
## 过滤在所有重复样本中小于1的基因，表达量太低也没研究意义
#count <- count[rowMeans(count)>1,]
##载入样本信息
fzdata2 <- read.csv('il17a注释信息.csv', sep = ',',header = T,row.names = 1)

#fzdata2 <- read.table("C:/Users/a/Desktop/表达矩阵/注释信息.txt",header = T,row.names = 1)

#一定要变为因子数据，否者用DESeq2包分析时候会出错
#fzdata2[,1] <- as.factor(fzdata2$type)

fazdata1 <- subset(fzdata2,type==1)
fazdata2 <- subset(fzdata2,type==2)
fazdata12 <- rbind(fazdata1,fazdata2)

as.data.frame(fazdata12)

all(rownames(fazdata12) %in% colnames(count))
all(rownames(fazdata12) == colnames(count))

fazdata12$type <- as.character(fazdata12$type)
fazdata12$patient_id <- as.character(fazdata12$patient_id)
dds <-  DESeqDataSetFromMatrix(countData = count,colData = fazdata12,design = ~type+patient_id)
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


write.csv(diff,kname1)

#Padj是P值矫正之后的数值，一般选取小于等于0.05（显著差异）的基因；同时log2FC是基因表达量的差异倍数。例如log2FC为1，证明这个基因在两种不同处理中的表达量相差了一倍，通常以大于1或小于-1为标准，大于1的为上调表达，少于-1的为下调表达。
foldChange = 1
padj = 0.05
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
dim(diffsig)

write.csv(diffsig,kname12)

up_DEG <- subset(diffsig, padj < 0.05 & log2FoldChange > 1)
down_DEG <- subset(diffsig, padj < 0.05 & log2FoldChange < -1)

write.csv(up_DEG,kname12up)
write.csv(down_DEG,kname12down)







