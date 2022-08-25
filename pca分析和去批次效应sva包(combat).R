

setwd("C:\\Users\\a\\Desktop\\Rrun1")
######表达矩阵批量读入与合并######
#先读入两个tsv检查gene_id和gene_name在不同样本中是否仍是一致的:
test1 = data.table::fread("C:\\Users\\a\\Desktop\\Rrun1\\AS_Bulk\\mRNA\\MA07_NK\\abundance.tsv")
test2 = data.table::fread("C:\\Users\\a\\Desktop\\Rrun1\\AS_Bulk\\mRNA\\MA14_CD4\\abundance.tsv")
identical(test1$gene_id,test2$gene_id)#返回逻辑值TRUE，确认一致
identical(test1$gene_name,test2$gene_name)#返回逻辑值TRUE，确认一致

#批量读取所有tsv格式后缀的文件名：
exp_all <- dir("C:\\Users\\a\\Desktop\\Rrun1\\AS_Bulk\\mRNA\\",
               pattern = "*.abundance.tsv$",#*表示任意前缀，$表示固定后缀
               recursive = T)
head(exp_all)

##注意，data.table::fread方式读入文件并不能指定行名，所以才使用这种方法
x1 <- data.table::fread("AS_Bulk\\mRNA\\MA07_NK\\abundance.tsv",
                        select = c("target_id"))
head(x1)




###这一行得改成用tximport的
5

#创建自定义函数,用于批量读入所有tcv文件的unstranded列（counts）：
exp_dt <- function(x){
  result <- data.table::fread(file.path("C:\\Users\\a\\Desktop\\Rrun1\\AS_Bulk\\mRNA\\",x),
                              select = c("unstranded"))
  return(result)
}
exp <- lapply(exp_all,exp_dt)#读取所有.tsv的unstranded列（将exp_all转换为list，并对list中的每一个元素都应用函数exp_dt）
exp <- do.call(cbind,exp)#将exp按列合并，并将list转化为data.table
exp[1:6,1:6]


#添加行名(gene_id),所有表达矩阵的gene_id是相同的：
exp <- as.data.frame(exp)#data.table格式不能自定行名，因此我们先转换为数据框
rownames(exp) <- x1$gene_id
exp[1:8,1:4]







 install.packages("tximport")
aa
 BiocManager::install("DESeq2")
 ?repositories