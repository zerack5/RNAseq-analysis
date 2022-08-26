setwd("C:\\Users\\a\\Desktop\\Rrun1")
library(tximport)
library(readr)
library(biomaRt)
library(rhdf5)

#creat a t2g file
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")


sample<- dir(file.path("C:/Users/a/Desktop/Rrun1/Control_Bulk/mRNA"))
files<-file.path("C:/Users/a/Desktop/Rrun1/Control_Bulk/mRNA",sample,"abundance.h5")
names(files)<-sample

#1 reading h5 file
#txi.kallisto<-tximport(files,type="kallisto",txOut=TRUE)
#txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")

#write.csv(txi.sum,"exp.csv",sep = "\t",quote = F,row.names = T,col.names = T)

#create t2g file --latest or read csv file downloaded
#library(biomaRt)
#mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#dataset = "hsapiens_gene_ensembl",
#host = 'ensembl.org')
#t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
# "external_gene_name"), mart = mart)
#t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
#ens_gene = ensembl_gene_id, ext_gene = external_gene_name)


tx2gene <- read.table("genename")

#2 reading tsv file
txi<-tximport(files,type = "kallisto",tx2gene = tx2gene,ignoreAfterBar = TRUE)
#txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")
write.csv(txi$counts,"exp2control.csv",sep = "\t",quote = F,row.names = T,col.names = T)


#222222222222222222222222222222222222222222222222222222222222222222
#creat a t2g file
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")


sample1<- dir(file.path("C:/Users/a/Desktop/Rrun1/AS_Bulk/mRNA"))
files1<-file.path("C:/Users/a/Desktop/Rrun1/AS_Bulk/mRNA",sample1,"abundance.h5")
names(files1)<-sample1

#1 reading h5 file
#txi.kallisto<-tximport(files,type="kallisto",txOut=TRUE)
#txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")

#write.csv(txi.sum,"exp.csv",sep = "\t",quote = F,row.names = T,col.names = T)

#create t2g file --latest or read csv file downloaded
#library(biomaRt)
#mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#dataset = "hsapiens_gene_ensembl",
#host = 'ensembl.org')
#t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
# "external_gene_name"), mart = mart)
#t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
#ens_gene = ensembl_gene_id, ext_gene = external_gene_name)


tx2gene <- read.table("genename")

#2 reading tsv file
txi1<-tximport(files1,type = "kallisto",tx2gene = tx2gene,ignoreAfterBar = TRUE)
#txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")
write.csv(txi1$counts,"exp2as.csv",sep = "\t",quote = F,row.names = T,col.names = T)


