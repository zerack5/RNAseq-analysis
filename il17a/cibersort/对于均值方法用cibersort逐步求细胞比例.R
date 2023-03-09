
twnamet1 <- "w1phase1.txt"
twnamet2 <- "w2phase2.txt"
twnamet3 <- "w3phase3.txt"


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

cresult1 <-"cibersortresult1.csv"
cresult2 <-"cibersortresult2.csv"
cresult3 <-"cibersortresult3.csv"


setwd(path0)


library(org.Hs.eg.db)

library("preprocessCore")

source('Cibersort.R')
result1 <- CIBERSORT("LM22.txt",twnamet1,perm = 1000, QN = T)  #perm 
write.csv(result1,file = cresult1)
 
source('Cibersort.R')
result2 <- CIBERSORT("LM22.txt",twnamet2,perm = 1000, QN = T)  #perm 
write.csv(result2,file = cresult2)


source('Cibersort.R')
result3 <- CIBERSORT("LM22.txt",twnamet3,perm = 1000, QN = T)  #perm 

write.csv(result3,file = cresult3)



