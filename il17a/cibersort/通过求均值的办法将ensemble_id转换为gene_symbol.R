




path0 <- "D:/数据分析中间文件/il17a/cibersort"
path1 <- "D:/数据分析中间文件/il17a/mRNAphase1"
path2 <- "D:/数据分析中间文件/il17a/mRNAphase2"
path3 <- "D:/数据分析中间文件/il17a/mRNAphase3"
sname1 <- "t1phase_il17a.txt"
sname1sum <- "t1phase_il17asum.txt"  
sname1sum <- "t1phase_il17asum.txt"
sname2 <- "t2phase_il17a.txt"
sname2sum <- "t2phase_il17asum.txt"
sname3 <- "t3phase_il17a.txt"
sname3sum <- "t3phase_il17asum.txt"


wnamet1 <- "w1phase_il17a.csv"
wnamet2 <- "w2phase_il17a.csv"
wnamet3 <- "w3phase_il17a.csv"

#install.packages("preprocessCore")
#BiocManager::install("preprocessCore")
setwd(path0)

#create t2g file --latest or read csv file downloaded
library(biomaRt)
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)


il17a1 <- read.table(sname1)

il17a1$ens_gene <- rownames(il17a1)


t3g <- unique(t2g[,-1])

il17a3 <- merge(il17a1,t3g,by="ens_gene")

il17a3_tmp <- il17a3
#il17a3_tmp$ext_gene <- as.factor(il17a3_tmp$ext_gene)


expr<- il17a3_tmp[,-1]

expr_mean=aggregate(.~ext_gene,mean,data=expr)

expr_mean <- expr_mean[-1,]

write.csv(expr_mean,file = wnamet1)


###########sname2


il17a1 <- read.table(sname2)

il17a1$ens_gene <- rownames(il17a1)


t3g <- unique(t2g[,-1])

il17a3 <- merge(il17a1,t3g,by="ens_gene")

il17a3_tmp <- il17a3
#il17a3_tmp$ext_gene <- as.factor(il17a3_tmp$ext_gene)


expr<- il17a3_tmp[,-1]

expr_mean=aggregate(.~ext_gene,mean,data=expr)

expr_mean <- expr_mean[-1,]

write.csv(expr_mean,file = wnamet2)



###########sname3


il17a1 <- read.table(sname3)

il17a1$ens_gene <- rownames(il17a1)


t3g <- unique(t2g[,-1])

il17a3 <- merge(il17a1,t3g,by="ens_gene")

il17a3_tmp <- il17a3
#il17a3_tmp$ext_gene <- as.factor(il17a3_tmp$ext_gene)


expr<- il17a3_tmp[,-1]

expr_mean=aggregate(.~ext_gene,mean,data=expr)

expr_mean <- expr_mean[-1,]

write.csv(expr_mean,file = wnamet3)





































#计算行平均值，按降序排列
index=order(rowMeans(expr[,-16]),decreasing = T)
#调整表达谱的基因顺序
expr_ordered=expr[index,]
#对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
keep=!duplicated(expr_ordered$genes)
#得到最后处理之后的表达谱矩阵
expr_max=expr_ordered[keep,]
expr_max




expr_max=aggregate(.~ext_gene,max,data=il17a3_tmp)



counts <- aggregate(il17a3_tmp, by = list(il17a3_tmp$X06_1_003), FUN=sum)


rownames(counts) <- counts$Group.1
tmp <- as.data.frame(table(counts$Group.1))
head(tmp[which(tmp$Freq >1),])
counts <-counts[,-1]
colnames(counts)





write.csv(il17a3,file="test.csv")



il17a1 <- read.table(sname1)
il17a1$ens_gene <- rownames(il17a1)
il17a1$ext_gene<-t2g[which(rownames(il17a1) %in% t2g$ens_gene),]$ext_gene


il17a1_tmp <- il17a1[,-c(15,16)]
counts <- aggregate(il17a1_tmp, by = list(il17a1$ext_gene), FUN=sum)
rownames(counts) <- counts$Group.1
tmp <- as.data.frame(table(counts$Group.1))
head(tmp[which(tmp$Freq >1),])
counts <-counts[,-1]
colnames(counts)