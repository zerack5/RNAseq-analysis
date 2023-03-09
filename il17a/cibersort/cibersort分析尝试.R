


path0 <- "D:/数据分析中间文件/il17a/cibersort"
path1 <- "D:/数据分析中间文件/il17a/mRNAphase1"
path2 <- "D:/数据分析中间文件/il17a/mRNAphase2"
path3 <- "D:/数据分析中间文件/il17a/mRNAphase3"
sname1 <- "t1phase_il17a.txt"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            sname1sum <- "t1phase_il17asum.txt"
sname2 <- "t2phase_il17a.txt"
sname2sum <- "t2phase_il17asum.txt"
sname3 <- "t3phase_il17a.txt"
sname3sum <- "t3phase_il17asum.txt"
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

tx2gene <- t2g

#BiocManager::install("wgcna")

library("wgcna")

il17b <- collapseRow()


lm22 <- read.table('LM22.txt', sep='\t')

il17a <- read.table(sname1,row.names = F)
rownames(il17a)<-t2g[which(rownames(il17a) %in% t2g$ens_gene),]$ext_gene



write.table(lm22,"LM22.txt")
lm22 <- read.csv("LM22.csv")
write.table(lm22,"LM22.txt")

library(org.Hs.eg.db)

library("preprocessCore")

source('Cibersort.R')
result1 <- CIBERSORT("LM22.txt",sname1,perm = 1000, QN = T)  #perm 
