getwd()
setwd("E:/drx/1/data4")

#换源
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
#


install.packages('Seurat')
install.packages("tidyverse")
install.packages("dplyr")
install.packaages("patchwork")



#
install.packages("sp")
install.packages("leidenbase")
BiocManager::install("leidenbase")
install.packages("car")
install.packages("carData")
#
BiocManager::install("glmGamPoi")

library(sp)
library(SeuratObject)
library(ggpubr)
library(leidenbase)
library(car)
library(carData)
library(survminer)

library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library("DoubletFinder")
library(glmGamPoi)
library(MAST)
library(SingleR)
library(celldex)
#查内存剩余
gc()

#读取10X数据
?Read10X
scRNA.counts=Read10X(data.dir="E:/")
class(scRNA.counts)

#创建Seurat对象
?CreateSeuratObject
scRNA=CreateSeuratObject(scRNA.counts,
                         min.cells=3,project="sampl.1",
                         min.features=40)
view(scRNA)

##质控
#查看细胞数量
table(scRNA@meta.data$orig.ident)
#计算质控指标
#计算线粒体基因比例
?PercentageFeatureSet
x=scRNA[[]]
###  ^：正则表达匹配开头；MT-：线粒体基因
scrna[["percentage.mt"]]=PercentageFeatureSet(scrna,pattern="^MT-")
#计算红细胞比例,先确定红细胞基因集
HB.genes<- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_M<- match(HB.genes,row.names(scrna@assays$RNA))
HB.genes<- row.names(scrna@assays$RNA)[HB_M]#找出目前数据集中，检测到的红细胞基因
HB.genes<-HB.genes[!is.na(HB.genes)]
scrna[["hb.percent"]]<-PercentageFeatureSet(scrna,features = HB.genes)

#可视化
violin<-VlnPlot(scrna,features = c("nFeature_RNA","nCount_RNA","percentage.mt","hb.percent"),
                pt.size = 0.01,ncol = 4)
#几个指标间相关性
plot1<- FeatureScatter(scrna,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot2<- FeatureScatter(scrna,feature1 = "nCount_RNA",feature2 = "percentage.mt")
plot3<- FeatureScatter(scrna,feature1 = "nCount_RNA",feature2 = "hb.percent")
#合并几个图象
polplot<- CombinePlots(plots = list(plot1,plot2,plot3),nrow=1,legend = "none")
#质控
scrna<- subset(scrna,subset = nFeature_RNA>500&percentage.mt<20&hb.percent<1&nCount_RNA>1000)
#标准化
scrna<-NormalizeData(scrna,normalization.method = "LogNormalize",scale.factor = 10000)

#保存图片
ggsave
#找寻高变基因
scrna1<-FindVariableFeatures(scrna,selection.method = "vst",nfeatures = 2000)
#数据中心化处理,内存大可以对全部基因进行处理，这里我们内存小，只选择高变基因
scale.genes<- VariableFeatures(scrna1)
scrna1<- ScaleData(scrna1,features = scale.genes)
#pca降维
scrna1<-RunPCA(scrna1,features = VariableFeatures(object = scrna1))
#聚类
scrna2<-FindNeighbors(scrna1,dims = 1:20)
scrna2<-FindClusters(scrna2,resolution=1)
scrna2<-BuildClusterTree(scrna2)
PlotClusterTree(scrna2)
#tsne降维
scrna2<-RunTSNE(scrna2,dims = 1:20)
embscrna2<-Embeddings(scrna2,'tsne')
write.csv(embscrna2,"embscrna2.csv")
plotn<- DimPlot(scrna2,reduction = "tsne")
ggsave(filename = "tsnereduction.pdf",width = 16,height = 9,dpi = 30000,plot = plotn)
#umap方式降维
scrna2<-RunUMAP(scrna2,dims = 1:20)
embscrna<-Embeddings(scrna2,'umap')
write.csv(embscrna,"umapscrna2.csv")
plotn1<-DimPlot(scrna2,reduction = "umap")
#常规分析，单个样本
kdo<- NormalizeData(kdo)
kdo <- FindVariableFeatures(kdo)
kdo<-ScaleData(kdo)
kdo <- RunPCA(kdo)
kdo<-FindNeighbors(kdo,dims = 1:30,reduction = "pca")
kdo<-FindClusters(kdo,resolution = 1)
kdo<-RunUMAP(kdo,dims = 1:30,reduction = "pca")



#自定义函数
noranddedou <- function(kdo){kdo<- NormalizeData(kdo)
kdo <- FindVariableFeatures(kdo)
kdo<-ScaleData(kdo)
kdo <- RunPCA(kdo)
nEXP<-round(0.08*length(kdo@meta.data$orig.ident)*length(kdo@meta.data$orig.ident)/10000)
kdod <-doubletFinder(kdo,PCs = 1:20,nExp = nEXP,pK = 0.09,pN = 0.25,sct = FALSE)#没用sct就false,第一次跑不要加reuse.pann，会报错
kdod$doublet.class<- kdod[[paste0("DF.classifications_0.25_0.09_",nEXP)]]
kdod[[paste0("DF.classifications_0.25_0.09_",2074)]]<-NULL
pann<-grep(pattern = "^pANN",x=names(kdod@meta.data),value = TRUE)
kdod$pANN_0.25_0.09_2074<-kdod[[pann]]
kdod[[pann]] <- NULL
kdod <-subset(kdod,subset = doublet.class != "Doublet")
return(kdod)
}
dirs<-list.dirs(path = "E:/prac/",recursive = FALSE,full.names = TRUE)
dirs<-subset(dirs,subset = dirs !=dirs[1])
kdod_l<-list()
for (i in seq_along(dirs)) {
  x<-Read10X(data.dir = dirs[i])
  x<-CreateSeuratObject(x,min.cells=3,min.features=40,project = "sample.1")
  x[["percentage.mt"]]=PercentageFeatureSet(x,pattern="^MT-")
  x<- subset(x,subset = nFeature_RNA>300&percentage.mt<20&nCount_RNA>1000)
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x,nFeatures=3000)
  x<-ScaleData(x)
  x<-RunPCA(x)
  nEXP<-round(0.08*length(x@meta.data$orig.ident)*length(x@meta.data$orig.ident)/10000)
  x <-doubletFinder(x,PCs = 1:30,nExp = nEXP,pK = 0.09,pN = 0.25,sct = FALSE)#没用sct就false,第一次跑不要加reuse.pann，会报错
  #去除
  x$doublet.class<- x[[paste0("DF.classifications_0.25_0.09_",nEXP)]]
  x[[paste0("DF.classifications_0.25_0.09_",nEXP)]]<-NULL
  pann<-grep(pattern = "^pANN",x=names(x@meta.data),value = TRUE)
  x[[paste0("paNN_0.25_0.09_",nEXP)]]<-x[[pann]]
  x[[pann]] <- NULL
  x <-subset(x,subset = doublet.class != "Doublet")
  kdod_l[[i]]<-x
  names(kdod_l)[i]<-paste0("c",i)
}
#############################################################
#############################################################
#############################################################
#############################################################
##########################################################################################################################
#############################################################
