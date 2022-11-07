
setwd("C:/Users/a/Desktop/表达矩阵")
b3gnt2as <- read.csv('exp2as.csv', sep = ',',header = T,row.names = 1)
b3gnt2c <- read.csv('exp2control.csv', sep = ',',header = T,row.names = 1)
nkb3gnt2as <- b3gnt2as[c(grep("_N",colnames(b3gnt2as)))]

a2 <- nkb3gnt2as["TNF",]

nkb3gnt2c <- b3gnt2c[c(grep("_N",colnames(b3gnt2c)))]
b2 <- nkb3gnt2c["TNF",]





setwd("C:/Users/a/Desktop/表达矩阵")
b3gnt2as <- read.csv('exp2as.csv', sep = ',',header = T,row.names = 1)
b3gnt2c <- read.csv('exp2control.csv', sep = ',',header = T,row.names = 1)
nkb3gnt2as <- b3gnt2as[c(grep("_M",colnames(b3gnt2as)))]

a2 <- nkb3gnt2as["TNF",]

nkb3gnt2c <- b3gnt2c[c(grep("_M",colnames(b3gnt2c)))]
b2 <- nkb3gnt2c["TNF",]





setwd("C:/Users/a/Desktop/表达矩阵")
b3gnt2as <- read.csv('exp2as.csv', sep = ',',header = T,row.names = 1)
b3gnt2c <- read.csv('exp2control.csv', sep = ',',header = T,row.names = 1)
nkb3gnt2as <- b3gnt2as[c(grep("CD4",colnames(b3gnt2as)))]

a2 <- nkb3gnt2as["TUBB",]

nkb3gnt2c <- b3gnt2c[c(grep("CD4",colnames(b3gnt2c)))]
b2 <- nkb3gnt2c["TUBB",]



setwd("C:/Users/a/Desktop/表达矩阵")
b3gnt2as <- read.csv('exp2as.csv', sep = ',',header = T,row.names = 1)
b3gnt2c <- read.csv('exp2control.csv', sep = ',',header = T,row.names = 1)
nkb3gnt2as <- b3gnt2as[c(grep("CD8",colnames(b3gnt2as)))]

a2 <- nkb3gnt2as["TUBB",]

nkb3gnt2c <- b3gnt2c[c(grep("CD8",colnames(b3gnt2c)))]
b2 <- nkb3gnt2c["TUBB",]



boxplot(nkb3gnt2as["ACTB",],nkb3gnt2c["B3GNT2",],main="NK细胞表达量",xlab="种类",ylab="表达量",col = c(2,3),names =c("AS","Normal"))






