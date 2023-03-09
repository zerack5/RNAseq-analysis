

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


setwd(path0)
library(tximport)
library(readr)
library(biomaRt)
library(rhdf5)

#creat a t2g file
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#k <- keys(txdb, keytype = "TXNAME")
#tx2gene <- select(txdb, k, "GENEID", "TXNAME")


sample<- dir(file.path(path1))
files<-file.path(path1,sample,"abundance.h5")
names(files)<-sample

#1 reading h5 file
#txi.kallisto<-tximport(files,type="kallisto",txOut=TRUE)
#txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")

#write.csv(txi.sum,"exp.csv",sep = "\t",quote = F,row.names = T,col.names = T)

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

#tx2gene <- read.table("genename")


#2 reading tsv file 
txi<-tximport(files,type = "kallisto",tx2gene = tx2gene,ignoreAfterBar = TRUE,ignoreTxVersion = TRUE)

txi.kallisto<-tximport(files,type="kallisto",txOut=TRUE)
txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")

write.table(txi$abundance,sname1,sep = "\t",quote = F,row.names = T,col.names = T)
write.table(txi.sum$abundance,sname1sum,sep = "\t",quote = F,row.names = T,col.names = T)




setwd(path0)
library(tximport)
library(readr)
library(biomaRt)
library(rhdf5)



sample<- dir(file.path(path2))
files<-file.path(path2,sample,"abundance.h5")
names(files)<-sample

#1 reading h5 file
#txi.kallisto<-tximport(files,type="kallisto",txOut=TRUE)
#txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")

#write.csv(txi.sum,"exp.csv",sep = "\t",quote = F,row.names = T,col.names = T)

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


#2 reading tsv file
txi<-tximport(files,type = "kallisto",tx2gene = tx2gene,ignoreAfterBar = TRUE,ignoreTxVersion = TRUE)
txi.kallisto<-tximport(files,type="kallisto",txOut=TRUE)
txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")
write.table(txi$abundance,sname2,sep = "\t",quote = F,row.names = T,col.names = T)
write.table(txi.sum$abundance,sname2sum,sep = "\t",quote = F,row.names = T,col.names = T)

########phase3

setwd(path0)
library(tximport)
library(readr)
library(biomaRt)
library(rhdf5)




sample<- dir(file.path(path3))
files<-file.path(path3,sample,"abundance.h5")
names(files)<-sample

#1 reading h5 file
#txi.kallisto<-tximport(files,type="kallisto",txOut=TRUE)
#txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")

#write.csv(txi.sum,"exp.csv",sep = "\t",quote = F,row.names = T,col.names = T)

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



#2 reading tsv file
txi<-tximport(files,type = "kallisto",tx2gene = tx2gene,ignoreAfterBar = TRUE,ignoreTxVersion = TRUE)
txi.kallisto<-tximport(files,type="kallisto",txOut=TRUE)
txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")
write.table(txi$abundance,sname3,sep = "\t",quote = F,row.names = T,col.names = T)
write.table(txi.sum$abundance,sname3sum,sep = "\t",quote = F,row.names = T,col.names = T)




