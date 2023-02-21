BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

BiocManager::install("tximport")
install.packages("ape")
install.packages("readr")
install.packages("rdf5")
BiocManager::install("rdf5")
install.packages("BiocManager")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rhdf5")



path0 <- "D:/数据分析中间文件/il17a"
path1 <- "D:/数据分析中间文件/il17a/mRNAphase1"
path2 <- "D:/数据分析中间文件/il17a/mRNAphase2"
path3 <- "D:/数据分析中间文件/il17a/mRNAphase3"
sname1 <- "1phase_il17a.csv"
sname2 <- "2phase_il17a.csv"
sname3 <- "3phase_il17a.csv"


setwd(path0)
library(tximport)
library(readr)
library(biomaRt)
library(rhdf5)

#creat a t2g file
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")


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

tx2gene <- read.table("genename")

#2 reading tsv file
txi<-tximport(files,type = "kallisto",tx2gene = tx2gene,ignoreAfterBar = TRUE,ignoreTxVersion = TRUE)
#txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")
write.csv(txi$abundance,sname1,sep = "\t",quote = F,row.names = T,col.names = T)





setwd(path0)
library(tximport)
library(readr)
library(biomaRt)
library(rhdf5)

#creat a t2g file
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")


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

tx2gene <- read.table("genename")

#2 reading tsv file
txi<-tximport(files,type = "kallisto",tx2gene = tx2gene,ignoreAfterBar = TRUE,ignoreTxVersion = TRUE)
#txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")
write.csv(txi$abundance,sname2,sep = "\t",quote = F,row.names = T,col.names = T)


########phase3

setwd(path0)
library(tximport)
library(readr)
library(biomaRt)
library(rhdf5)

#creat a t2g file
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")


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

tx2gene <- read.table("genename")

#2 reading tsv file
txi<-tximport(files,type = "kallisto",tx2gene = tx2gene,ignoreAfterBar = TRUE,ignoreTxVersion = TRUE)
#txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")
write.csv(txi$abundance,sname3,sep = "\t",quote = F,row.names = T,col.names = T)




