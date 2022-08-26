

dataas <- read.csv('expas.csv', sep = ',',header = T,row.names = 1)
datac <- read.csv('expcontrol.csv', sep = ',',header = T,row.names = 1)
cbdata <- cbind(dataas, datac)


fzdata <- read.csv('注释信息.csv', sep = ',',header = T,row.names = 1)

model <- model.matrix(~as.factor(fzdata$type))


############limma矫正
#modcombat <- model.matrix(as.formula(paste('~', fzdata_type, sep=" ")), data=fzdata)
library(limma)
suppressMessages(library(DESeq2))
suppressMessages(library("RColorBrewer"))
)
suppressMessages(library("amap"))
suppressMessages(library("ggplot2"))
suppressMessages(library("BiocParallel"))
suppressMessages(library("YSX"))
suppressMessages(library(sva))
suppressMessages(library(ggfortify))
suppressMessages(library(patchwork))
suppressMessages(library(ggbeeswarm))
suppressMessages(library(limma))

#BiocManager::install("amap")
#BiocManager::install("BiocParallel")
#BiocManager::install("patchwork")
#BiocManager::install("ggbeeswarm")
BiocManager::install("YSX")
BiocManager::install("gplots")

modcombat <- model

library("ggplot2")
SV = fzdata$batch
expr_mat_batch_correct_limma1 <- removeBatchEffect(cbdata, covariates = SV, design=modcombat)


setwd("C:/Users/a/Desktop/表达矩阵")
middleoption2 <- data.frame(expr_mat_batch_correct_limma1)

select2 <- middleoption2[c(grep("_M",colnames(middleoption2)))]
f <- select2["GAPDH",]


library("devtools")
#install_github("kassambara/factoextra")
library("FactoMineR")
library("factoextra")
af1.pca <- PCA(t(expr_mat_batch_correct_limma1),graph = FALSE)
combat.pca <- PCA(t(expr_mat_batch_correct_limma1),graph = FALSE)
fviz_pca_ind(combat.pca,
             geom= "point",
             col.ind = fzdata$label,
             addEllipses = TRUE,
             legend.title="Group"  )




#sp_pca(expr_mat_batch_correct_limma1[1:5000,], cbdata,
#color_variable="fzdata$type", shape_variable = "fzdata$label") +
#  aes(size=1) + guides(size = F)
#?sp_pca

#expr_mat_batch_correct_limma1


########################查找矫正之前的值
setwd("C:/Users/a/Desktop/表达矩阵")
selectascb <- cbdata[c(grep("_M",colnames(cbdata)))]
d <- selectascb["GAPDH",]



