




dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design= ~ batch + condition) #~在R里面用于构建公式对象，~左边为因变量，右边为自变量。
dds <- DESeq(dds) #标准化
res <- results(dds, contrast=c("condition","treated","control")) #差异分析结果
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design= ~ batch + condition) #~在R里面用于构建公式对象，~左边为因变量，右边为自变量。


#过滤低丰度数据 dds <- dds[rowSums(counts(dds)) > 1, ] 

vsd <- vst(object=dds,blind=T) 

rld <- rlog(object=dds,blind=F) 



