#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

setwd("C:\\Users\\a\\Desktop\\Rrun1")
library(tximport)
library(readr)
library(biomaRt)
library(rhdf5)



library(tximportData)
dir <- system.file("extdata", package = "tximportData")
list.files(dir)


#creat a t2g file
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#k <- keys(txdb, keytype = "TXNAME")
#tx2gene <- select(txdb, k, "GENEID", "TXNAME")


library(readr)
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
head(tx2gene)
  

library(tximport)

#txi<-tximport(files,type = "salmon", tx2gene = tx2gene)
#names(txi)
  
  
sample<- dir(file.path("C:/Users/a/Desktop/Rrun1/AS_Bulk/mRNA"))
files<-file.path("C:/Users/a/Desktop/Rrun1/AS_Bulk/mRNA",sample,"abundance.h5")
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


#tx2gene <- read.table("genename")

#2 reading tsv file
txi<-tximport(files,type = "kallisto",tx2gene = tx2gene,ignoreAfterBar = TRUE)
#txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")
write.csv(txi$abundance,"expas.csv",sep = "\t",quote = F,row.names = T,col.names = T)


