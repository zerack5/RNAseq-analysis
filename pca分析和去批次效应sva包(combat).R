#####################

## 4. DESeq2包消除批次效应
# design 设计矩阵中加入引起批次效应的因素(SR or PE)

cbmdata <-as.matrix(cbdata)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cbmdata,
                              colData = colData,
                              design = ~ condition+ type)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="condition_untreated_vs_treated")

summary(res)
head(res)
resOdered <- res[order(res$padj),]
deg <- as.data.frame(resOdered)
#deg <- na.omit(deg)
dim(deg)
#write.csv(deg,file= "diff_deseq2.csv")

# 看看批次效应导致的差异表达基因,
# 这里是测序方法SR(single-end)和PE (paired-end)
batch_res <-  results(dds, name="type_SR_vs_PE")
summary(batch_res)
head(batch_res)
batchResOdered <- res[order(batch_res$padj),]
batchDeg <- as.data.frame(batchResOdered)
batchDeg <- na.omit(batchDeg)
dim(batchDeg)
head(batchDeg)


###################################################












#BiocManager::install("sva")
#BiocManager::install("bladderbatch")

#建立批次效应的模型，data$type表示的是数据中除了有不同的批次，还有生物学上的差异。
#在校正的时候要保留生物学上的差异，不能矫正过枉。
dataas <- read.csv('expas.csv', sep = ',',header = T,row.names = 1)
datac <- read.csv('expcontrol.csv', sep = ',',header = T,row.names = 1)
cbdata <- cbind(dataas, datac)


fzdata <- read.csv('注释信息.csv', sep = ',',header = T,row.names = 1)

model <- model.matrix(~as.factor(fzdata$type))

library("devtools")
#install_github("kassambara/factoextra")
library("FactoMineR")
library("factoextra")
af1.pca <- PCA(t(cbdata),graph = FALSE)
combat.pca <- PCA(t(cbdata),graph = FALSE)
fviz_pca_ind(combat.pca,
             geom= "point",
             col.ind = fzdata$label,
             addEllipses = TRUE,
             legend.title="Group"  )

##消除批次效应
library(sva)

combat_Expr <- ComBat(dat = cbdata,batch = fzdata$batch,mod = model)

##combat_Expr就是校正后的数据
#BiocManager::install("FactoMineR")
 
#BiocManager::install("factoextra")
#install.packages("devtools")

library("devtools")
#install_github("kassambara/factoextra")
library("FactoMineR")
library("factoextra")
af1.pca <- PCA(t(combat_Expr),graph = FALSE)
combat.pca <- PCA(t(combat_Expr),graph = FALSE)
fviz_pca_ind(combat.pca,
             geom= "point",
             col.ind = fzdata$label,
             addEllipses = TRUE,
             legend.title="Group"  )
