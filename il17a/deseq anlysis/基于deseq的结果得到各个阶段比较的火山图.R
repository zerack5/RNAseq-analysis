
#install.packages("pheatmap")

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')


if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('ggrepel')


cphase12 <- "The difference between phase1 and phase2"
cphase23 <- "The difference between phase2 and phase3"



#########volcano plot of deseq result of phase1 and phase2


heat1data<-read.csv("D:/数据分析中间文件/il17a/all_diff12_paired_gene_res.csv", header = TRUE)


colnames(heat1data)
heatdata <- heat1data[,c("X","log2FoldChange","pvalue","padj" )]



library(EnhancedVolcano)
EnhancedVolcano(heat1data,
                lab = heat1data$X,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'difference between phase1 and phase2',
                pCutoff = 10e-4,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0)






#heatdata <- heatdata[which(abs(heatdata$log2FoldChange)>1),]
library(ggrepel)
library(ggplot2)
ggplot(heatdata,aes(log2FoldChange,-log10(padj)))+ 
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)")) #修改坐标轴命名细节



data <- heatdata
data$changed <- factor(ifelse(data$padj < 0.001 & abs(data$log2FoldChange) > 1, 
                              ifelse(data$log2FoldChange > 1,'positive','negative'),'no'))
data$changed <- as.character(data$changed)

data$selectedgene <- ifelse(data$padj<0.001 & abs(data$log2FoldChange)>1,data$X,NA)
head(data)

#概览一下选出的gene数量：
summary(data$changed)

ggplot(data=data,aes(x =log2FoldChange,y = -log10(padj),
                #分别给正负显著变化的基因在图中根据颜色、大小标注出来：
                color=factor(changed),
                size=factor(changed)))+  
  geom_point(alpha=0.8 , size =1 ,shape =19)+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)"))+
  theme_grey(base_size = 15)+
  ggtitle(cphase12)+
  
  scale_color_manual(values = c('blue','grey','red'))+
  scale_size_manual(values = c(2,1,2))+
  theme(panel.background = element_rect(fill = "gray90", color = "white"),#改变图表面板（panel）的背景颜色和样式
        panel.grid.major = element_blank(),#s删除散点图的主要网格线
        panel.grid.minor = element_blank(),#删除散点图中的次要网格线
        legend.title = element_blank(), #图例的设置参数
        legend.position = c(0.1,0.9),
        legend.background = element_rect(fill='transparent'))+

geom_hline(yintercept=-log10(0.001),linetype=4)+
geom_vline(xintercept=c(-1,1),linetype=4)+


geom_text_repel(aes(label=selectedgene), color="black",size=3,#为显著的基因基于算法不重合地添加genename
                    box.padding=unit(0.5, "lines"), 
                    point.padding=NA, 
                    segment.colour = "black") 

ggsave("The difference between phase1 and phase2.tiff",device = "tiff",scale = 2,width = NA,height = NA,dpi = 300,path = "D:/数据分析中间文件/il17a/figure")




#########volcano plot of deseq result of phase2 and phase3



heat2data<-read.csv("D:/数据分析中间文件/il17a/all_diffsig23_paired_gene_res.csv", header = TRUE)


colnames(heat2data)
heatdata <- heat2data[,c("X","log2FoldChange","pvalue","padj" )]



library(EnhancedVolcano)
EnhancedVolcano(heat2data,
                lab = heat2data$X,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'difference between phase1 and phase2',
                pCutoff = 10e-4,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0)






#heatdata <- heatdata[which(abs(heatdata$log2FoldChange)>1),]
library(ggrepel)
library(ggplot2)
ggplot(heatdata,aes(log2FoldChange,-log10(padj)))+ 
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)")) #修改坐标轴命名细节



data <- heatdata
data$changed <- factor(ifelse(data$padj < 0.001 & abs(data$log2FoldChange) > 1, 
                              ifelse(data$log2FoldChange > 1,'positive','negative'),'no'))
data$changed <- as.character(data$changed)

data$selectedgene <- ifelse(data$padj<0.001 & abs(data$log2FoldChange)>1,data$X,NA)
head(data)

#概览一下选出的gene数量：
summary(data$changed)

ggplot(data=data,aes(x =log2FoldChange,y = -log10(padj),
                     #分别给正负显著变化的基因在图中根据颜色、大小标注出来：
                     color=factor(changed),
                     size=factor(changed)))+  
  geom_point(alpha=0.8 , size =1 ,shape =19)+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)"))+
  theme_grey(base_size = 15)+
  ggtitle(cphase23)+
  
  scale_color_manual(values = c('blue','grey','red'))+
  scale_size_manual(values = c(2,1,2))+
  theme(panel.background = element_rect(fill = "gray90", color = "white"),#改变图表面板（panel）的背景颜色和样式
        panel.grid.major = element_blank(),#s删除散点图的主要网格线
        panel.grid.minor = element_blank(),#删除散点图中的次要网格线
        legend.title = element_blank(), #图例的设置参数
        legend.position = c(0.1,0.9),
        legend.background = element_rect(fill='transparent'))+
  
  geom_hline(yintercept=-log10(0.001),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  
  
  geom_text_repel(aes(label=selectedgene), color="black",size=3,#为显著的基因基于算法不重合地添加genename
                  box.padding=unit(0.5, "lines"), 
                  point.padding=NA, 
                  segment.colour = "black") 



ggsave("the difference between phase2 and phase3.tiff",device = "tiff",scale = 2,width = NA,height = NA,dpi = 300,path = "D:/数据分析中间文件/il17a/figure")


