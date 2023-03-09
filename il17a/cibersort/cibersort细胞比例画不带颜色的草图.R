
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
cresult2 <-"cibersortresult3.csv"






setwd(path0)
#读取数据sampleas
sample1 <- read.csv(cresult1, sep = ',',header = T,row.names = 1)
sample2 <- read.csv(cresult2, sep = ',',header = T,row.names = 1)
sample3 <- read.csv(cresult3, sep = ',',header = T,row.names = 1)



#subset提取相应的NK行列
nk1 <- subset(sample1,MARK=="NK",select = c("NK.cells.resting","NK.cells.activated","Monocytes"))

nk2 <- subset(sample2,MARK=="NK",select = c("NK.cells.resting","NK.cells.activated","Monocytes"))

nk3 <- subset(sample3,MARK=="NK",select = c("NK.cells.resting","NK.cells.activated","Monocytes"))



#subset提取相应的Mo行列
Mo1 <- subset(sample1,MARK=="Mo",select = c("Monocytes","Macrophages.M0","Macrophages.M1","Macrophages.M2"))
Mo2 <- subset(sample2,MARK=="Mo",select = c("Monocytes","Macrophages.M0","Macrophages.M1","Macrophages.M2"))
Mo3 <- subset(sample3,MARK=="Mo",select = c("Monocytes","Macrophages.M0","Macrophages.M1","Macrophages.M2"))


Moas
Moc

#subset提取相应的CD4行列  
CD41 <- subset(sample1,MARK=="CD4",select = c("T.cells.CD8","T.cells.CD4.naive","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper","T.cells.regulatory..Tregs.","T.cells.gamma.delta"))
CD42 <- subset(sample2,MARK=="CD4",select = c("T.cells.CD8","T.cells.CD4.naive","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper","T.cells.regulatory..Tregs.","T.cells.gamma.delta"))
CD43 <- subset(sample3,MARK=="CD4",select = c("T.cells.CD8","T.cells.CD4.naive","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper","T.cells.regulatory..Tregs.","T.cells.gamma.delta"))



#subset提取相应的CD8行列
CD81 <- subset(sample1,MARK=="CD8",select = c("T.cells.CD8","T.cells.CD4.naive","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper","T.cells.regulatory..Tregs.","T.cells.gamma.delta"))
CD82 <- subset(sample2,MARK=="CD8",select = c("T.cells.CD8","T.cells.CD4.naive","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper","T.cells.regulatory..Tregs.","T.cells.gamma.delta"))
CD83 <- subset(sample3,MARK=="CD8",select = c("T.cells.CD8","T.cells.CD4.naive","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper","T.cells.regulatory..Tregs.","T.cells.gamma.delta"))


#ggplot图NK
boxplot(sample1$NK.cells.resting,sample2$NK.cells.resting,sample3$NK.cells.resting,main="NK细胞比例",xlab="细胞种类",ylab="比例")

boxplot(sample1$NK.cells.activated,sample2$NK.cells.activated,sample3$NK.cells.activated,main="NK细胞比例",xlab="细胞种类",ylab="比例")

boxplot(sample1$B.cells.naive,sample2$B.cells.naive,sample3$B.cells.naive,main="B.naive细胞比例",xlab="细胞种类",ylab="比例")


sample1



#ggplot图Mo
boxplot(sample1$Monocytes,sample2$Monocytes,sample3$Monocytes,main="Monocytes细胞比例",xlab="细胞种类",ylab="比例")
boxplot(sample1$Macrophages.M0,sample2$Macrophages.M0,sample3$Macrophages.M0,main="Macrophages.M0细胞比例",xlab="细胞种类",ylab="比例")
boxplot(sample1$Macrophages.M1,sample2$Macrophages.M1,sample3$Macrophages.M1,main="Macrophages.M1细胞比例",xlab="细胞种类",ylab="比例")
boxplot(sample1$Macrophages.M2,sample2$Macrophages.M2,sample3$Macrophages.M2,main="Macrophages.M2细胞比例",xlab="细胞种类",ylab="比例")





#ggplot图Bcells
boxplot(sample1$B.cells.memory,sample2$B.cells.memory,sample3$B.cells.memory,main="B.cells.memory细胞比例",xlab="细胞种类",ylab="比例")
boxplot(sample1$Plasma.cells,sample2$Plasma.cells,sample3$Plasma.cells,main="Plasma.cells比例",xlab="细胞种类",ylab="比例")


#ggplot图Tcells
boxplot(sample1$T.cells.CD8,sample2$T.cells.CD8,sample3$T.cells.CD8,main="T.cells.CD8比例",xlab="细胞种类",ylab="比例")
boxplot(sample1$T.cells.CD4.naive,sample2$T.cells.CD4.naive,sample3$T.cells.CD4.naive,main="T.cells.CD4.naive比例",xlab="细胞种类",ylab="比例")

boxplot(sample1$T.cells.CD4.memory.resting,sample2$T.cells.CD4.memory.resting,sample3$T.cells.CD4.memory.resting,main="T.cells.CD4.memory.resting比例",xlab="细胞种类",ylab="比例")

boxplot(sample1$T.cells.CD4.memory.activated,sample2$T.cells.CD4.memory.activated,sample3$T.cells.CD4.memory.activated,main="T.cells.CD4.memory.activated比例",xlab="细胞种类",ylab="比例")

boxplot(sample1$T.cells.follicular.helper,sample2$T.cells.follicular.helper,sample3$T.cells.follicular.helper,main="T.cells.follicular.helper比例",xlab="细胞种类",ylab="比例")

boxplot(sample1$T.cells.regulatory..Tregs.,sample2$T.cells.regulatory..Tregs.,sample3$T.cells.regulatory..Tregs.,main="T.cells.regulatory..Tregs.比例",xlab="细胞种类",ylab="比例")

boxplot(sample1$T.cells.gamma.delta,sample2$T.cells.gamma.delta,sample3$T.cells.gamma.delta,main="T.cells.gamma.delta比例",xlab="细胞种类",ylab="比例")



boxplot(sample1$Dendritic.cells.resting,sample2$Dendritic.cells.resting,sample3$Dendritic.cells.resting,main="Dendritic.cells.resting比例",xlab="细胞种类",ylab="比例")

boxplot(sample1$Dendritic.cells.activated,sample2$Dendritic.cells.activated,sample3$Dendritic.cells.activated,main="Dendritic.cells.activated比例",xlab="细胞种类",ylab="比例")


boxplot(sample1$Mast.cells.resting,sample2$Mast.cells.resting,sample3$Mast.cells.resting,main="Mast.cells.resting比例",xlab="细胞种类",ylab="比例")


boxplot(sample1$Mast.cells.activated,sample2$Mast.cells.activated,sample3$Mast.cells.activated,main="Mast.cells.activated比例",xlab="细胞种类",ylab="比例")



boxplot(sample1$Eosinophils,sample2$Eosinophils,sample3$Eosinophils,main="Eosinophils比例",xlab="细胞种类",ylab="比例")


boxplot(sample1$Neutrophils,sample2$Neutrophils,sample3$Neutrophils,main="Neutrophils比例",xlab="细胞种类",ylab="比例")


###############################################################################################################################






boxplot(CD4as$T.cells.CD8,CD4c$T.cells.CD8,main="CD4细胞比例",xlab="组别",ylab="比例")
boxplot(CD4as$T.cells.CD4.naive,CD4c$T.cells.CD4.naive,main="CD4细胞比例",xlab="组别",ylab="比例")
boxplot(CD4as$T.cells.CD4.memory.activated,CD4c$T.cells.CD4.memory.activated,main="CD4细胞比例",xlab="组别",ylab="比例")
boxplot(CD4as$T.cells.CD4.memory.resting,CD4c$T.cells.CD4.memory.resting,main="CD4细胞比例",xlab="组别",ylab="比例")
boxplot(CD4as$T.cells.follicular.helper,CD4c$T.cells.follicular.helper,main="CD4细胞比例",xlab="组别",ylab="比例")
boxplot(CD4as$T.cells.regulatory..Tregs,CD4c$T.cells.regulatory..Tregs,main="CD4细胞比例",xlab="组别",ylab="比例")
boxplot(CD4as$T.cells.gamma.delta,CD4c$T.cells.gamma.delta,main="CD4细胞比例",xlab="组别",ylab="比例")


#ggplot图CD8
boxplot(CD8as$T.cells.CD8,CD4c$T.cells.CD8,main="CD8细胞比例",xlab="组别",ylab="比例")
boxplot(CD8as$T.cells.CD4.naive,CD4c$T.cells.CD4.naive,main="CD8细胞比例",xlab="组别",ylab="比例")
boxplot(CD8as$T.cells.CD4.memory.activated,CD4c$T.cells.CD4.memory.activated,main="CD8细胞比例",xlab="组别",ylab="比例")
boxplot(CD8as$T.cells.CD4.memory.resting,CD4c$T.cells.CD4.memory.resting,main="CD8细胞比例",xlab="组别",ylab="比例")
boxplot(CD8as$T.cells.follicular.helper,CD4c$T.cells.follicular.helper,main="CD8细胞比例",xlab="组别",ylab="比例")
boxplot(CD8as$T.cells.regulatory..Tregs,CD4c$T.cells.regulatory..Tregs,main="CD8细胞比例",xlab="组别",ylab="比例")
boxplot(CD8as$T.cells.gamma.delta,CD4c$T.cells.gamma.delta,main="CD8细胞比例",xlab="组别",ylab="比例")



















